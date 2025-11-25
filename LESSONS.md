# Lessons Learned: Building Torsor in One Night

**Date:** Saturday, November 22, 2025
**Location:** Canada (Saturday night hacking session)
**Goal:** Build a parametric CAD system inspired by MathCAD that combines algebraic geometry with analytical stress analysis

---

## The Initial Vision

**The Request:** "My boss wants to build a better CAD product. I have C++ expertise from Google HEIR, background in sheaf theory, and experience with industrial automation (Orin/Lexium drives). Can we build something like MathCAD but for 3D?"

**The Constraint:** "STL files" kept coming up as a key technology.

**The Insight:** STL isn't the *source format* - it's the *compilation target*. Just like how you write C++ source code and compile to machine code, Torsor writes algebraic geometry and "compiles" to STL for consumption by other tools.

---

## Design Evolution: From Meshes to Algebra

### The Traditional CAD Approach (B-Rep/NURBS)
- Boundary Representation: Store surfaces as triangle meshes
- Parametric changes require **remeshing** (expensive)
- Stress analysis requires **FEA** (discretize, solve linear system)
- Iteration loop: Design → Mesh → Simulate → Redesign

### The Torsor Approach (F-Rep/SDF)
- Function Representation: Geometry is an **algebraic function** f(x,y,z)
- Signed Distance Field: `f(p) < 0` means "inside", `f(p) > 0` means "outside"
- Parametric changes are **instant** (just function parameters)
- Stress analysis uses **Roark's Formulas** (analytical, no meshing for simple cases)
- No iteration: Design → Calculate → Done

### Why This Works
**Sheaf Theory Mapping:**
- Local Sections = Geometric primitives (Sphere, Cylinder, Beam)
- Gluing Maps = Boolean operations (Union, Intersection, Difference)
- Stalks = Point evaluations f(x,y,z)

**The Killer Insight:** A joint isn't a separate entity (nut, bolt, weld) - it's the **geometric intersection** of two bodies. `Intersection(beam, support)` directly gives you the contact patch.

---

## Technical Decisions

### Decision 1: C++20 (Not C++23)
**Why:** Mac's Clang doesn't fully support C++23 (`-std=c++23` fails, requires `-std=c++2b`)
**Consequence:** Used C++20 concepts, still got all the template magic we needed
**Lesson:** Target the compiler version your users actually have

### Decision 2: Dual Numbers for Autodiff
**What:** `Dual<T> = {value, gradient}` enables automatic differentiation
**Why:** Forward-mode autodiff gives exact derivatives with minimal code
**Implementation:**
```cpp
Dual operator*(const Dual& rhs) const {
    // Product rule: (fg)' = f'g + fg'
    return {val * rhs.val, grad * rhs.val + val * rhs.grad};
}
```
**Limitation:** Only tracks one directional derivative at a time (not the full gradient vector)

### Decision 3: Concepts for Type Safety
**Pattern:**
```cpp
template<typename T>
concept Scalar = std::floating_point<T> || requires(T t) {
    { t.val } -> std::floating_point;
    { t.grad } -> std::floating_point;
};

template<typename G, typename T>
concept ImplicitField = Scalar<T> && requires(const G& g, std::array<T,3> p) {
    { g.eval(p) } -> std::convertible_to<T>;
};
```
**Result:** All primitives, boolean ops, and rendering work with both `double` and `Dual<double>` seamlessly

### Decision 4: Binary PPM Output
**Initial:** P3 ASCII format (4.1MB, incompatible with ParaView)
**Fix:** P6 binary format (1.4MB, ParaView reads it)
**Lesson:** Check your file formats against the tools users actually have

---

## Architecture Insights

### Backend Agnostic Design
**Separation of Concerns:**
1. **Geometry Layer** (platform-independent)
   - Primitives define `eval()` templated on Scalar type
   - Boolean ops compose via templates
   - No CPU/GPU specifics here

2. **Evaluation Engine** (platform-specific)
   - CPU version: Direct function calls
   - GPU version (future): CUDA/Metal kernels
   - WebGPU version (future): WGSL shaders

**Why This Matters:** Write geometry once, run anywhere (Mac dev → Orin prod)

### The Renderer vs. The Solver
This is the **critical performance lesson** from tonight.

#### For Rendering (Raymarching)
**Task:** Compute normals for lighting (millions per frame)
**Current Approach:** Finite differences
```cpp
fx = (f(x+ε) - f(x-ε)) / 2ε  // 6 evaluations total
```
**Why It's Right:**
- Passes simple `double` types (8 bytes, fits in registers)
- Instruction cache friendly (tight loop, no polymorphism)
- "Good enough" accuracy for visualization
- ARM64 loves this (no register spills)

**Wrong Approach:** Using `Dual3<double>` (32 bytes) would cause register pressure and stack spills

#### For Optimization (Future Solver)
**Task:** Minimize weight subject to stress < yield (once per iteration)
**Future Approach:** `Gradient<N>` for parameter derivatives
```cpp
Gradient<double, 2> r = {5.0, {1.0, 0.0}};  // Seed ∂/∂r
Gradient<double, 2> h = {10.0, {0.0, 1.0}}; // Seed ∂/∂h
auto stress = calculate_stress(r, h, load);
// stress.grad[0] = ∂σ/∂r EXACTLY (no epsilon error)
```
**Why It's Right:**
- Need exact gradients for gradient descent convergence
- Only compute once per design iteration (not millions of times)
- The 32-byte overhead is justified by algorithmic superiority

**Key Distinction:**
- **Spatial derivatives** (normals): ∇f(x,y,z) → Use finite differences
- **Parameter derivatives** (optimization): ∂σ/∂r → Use autodiff

---

## Engineering Victories

### Victory 1: Parametric Recalculation
Change `beam.depth = 6.0` to `beam.depth = 8.0`:
- Section properties update instantly (pure algebra)
- Stresses recalculate analytically (Roark's formulas)
- No remeshing, no FEA iteration

**This is MathCAD for 3D.**

### Victory 2: Analytical Stress Without FEA
**Traditional Approach:**
1. Model geometry in CAD
2. Export to FEA tool
3. Mesh the geometry (100k+ elements)
4. Solve Ku=F linear system (expensive)
5. Extract stresses from mesh

**Torsor Approach:**
1. Define geometry algebraically
2. Calculate section properties (A, I, S) from geometry
3. Apply Roark's formula: σ = Mc/I
4. Done (instant)

**When This Works:** Simple structural loading (beams, contact stress)
**When You Need FEA:** Complex geometries, nonlinear materials, dynamic loads

### Victory 3: C6x13 Beam Demo
**Results:**
- Beam: C6x13 channel (24" long) on 2" diameter cylinder
- Load: 1000 lbs point load at center
- Bending stress: 974.8 psi
- Contact stress: 5527 psi (Hertzian)
- Safety factor: 36.9 (A36 steel, yield = 36 ksi)
- **Status: DESIGN OK** ✓

**What This Proves:** The entire workflow works end-to-end.

---

## Errors and Fixes

### Error 1: C++23 Standard Not Available
**Error:** `invalid value 'c++23' in '-std=c++23'`
**Fix:** Changed to `-std=c++20`
**Root Cause:** Mac Clang uses `-std=c++2b` for experimental C++23

### Error 2: min/max Ambiguity
**Error:** `no matching function for call to 'min'`
**Fix:** Added overloads for plain `double`:
```cpp
template <std::floating_point T>
constexpr T min(T a, T b) { return (a < b) ? a : b; }
```
**Root Cause:** Boolean ops needed to work with both `Dual<T>` and plain `T`

### Error 3: Union Type Deduction Failed
**Error:** `no viable constructor or deduction guide`
**Fix:** Explicit template args:
```cpp
Union<decltype(beam_geom), Cylinder> scene{beam_geom, support};
```
**Root Cause:** C++20 CTAD couldn't infer nested template types

### Error 4: PPM Format Wrong
**Error:** ParaView couldn't read P3 ASCII format
**Fix:** Switched to P6 binary format
**Lesson:** Always test with the actual tools users will use

---

## The Sheaf Theory Connection

**From the research papers shared:**
- Sheaf-Wreath Attention for algebraic learning
- Galois-Sheaf theory for categorical bootstrap
- Financial arbitrage via sheaf cohomology

**Applied to Torsor:**

### Local Sections (Primitives)
Each primitive (Sphere, Cylinder, Beam) is a **section** over 3D space:
```cpp
f: ℝ³ → ℝ  (the signed distance function)
```

### Gluing Maps (Boolean Operations)
Union, Intersection, Difference are **gluing constraints**:
```cpp
Union(f, g)(p) = min(f(p), g(p))
Intersection(f, g)(p) = max(f(p), g(p))
```

### The Stalks (Point Evaluations)
At each point p in space, the "stalk" is the value f(p) telling you distance to surface.

### Why This Matters
Sheaf theory guarantees **consistency** - if two primitives agree on their overlap, the glued result is well-defined. This is why boolean operations "just work" without special case handling.

---

## What We Shipped Tonight

**Files:**
- `main.cpp` (847 lines): Complete implementation
- `README.md`: Documentation and theory
- `.gitignore`: Clean repo management
- `LESSONS.md`: This file

**Core Features:**
- Dual numbers (forward-mode autodiff)
- Primitives: Sphere, Box, Cylinder, Beam
- Boolean ops: Union, Intersection, Difference
- C6x13 channel beam (AISC standard)
- Section property calculations
- Roark's formulas (bending, contact stress)
- Raymarching renderer
- Binary PPM output

**Engineering Results:**
- Working parametric geometry kernel
- Analytical stress calculations
- Safety factor analysis
- 3D visualization

---

## Future Roadmap

### Phase 1: The Optimizer (Near-term)
**Goal:** Minimize weight subject to stress < yield

**Implementation:**
```cpp
template <typename T, size_t N>
struct Gradient {
    T val;
    std::array<T, N> grad;  // Jacobian vector
    // ... implement arithmetic ...
};

// Optimize a cylinder (r, h are variables)
using Var = Gradient<double, 2>;
Var r = {5.0, {1.0, 0.0}};  // ∂/∂r
Var h = {10.0, {0.0, 1.0}}; // ∂/∂h

auto stress = calculate_stress(r, h, load);
// stress.grad[0] = ∂σ/∂r exactly
// Use gradient descent to minimize volume
```

**This Beats SolidWorks:**
- SolidWorks: Run FEA 100 times with perturbed parameters (slow, approximate)
- Torsor: Get exact gradients ∂σ/∂r in one pass (fast, exact)

### Phase 2: More Geometry (Medium-term)
**AISC Structural Shapes:**
- I-beams (W-shapes)
- Angles (L-shapes)
- Tees (WT-shapes)
- Hollow structural sections (HSS)

**More Roark's Formulas:**
- Chapter 9: Curved beams
- Chapter 10: Torsion
- Chapter 13: Shells of revolution
- Chapter 15: Columns

### Phase 3: GPU Acceleration (Medium-term)
**Mac (Development):**
- Metal compute shaders
- Raymarching on GPU (interactive framerates)

**NVIDIA Orin (Production):**
- CUDA kernels
- Real-time optimization on robot

### Phase 4: Advanced Features (Long-term)
**STL Import → Algebraic Fitting:**
- Read triangle mesh
- Fit SDF using neural network or optimization
- Get parametric representation of scanned parts

**WebGPU Version:**
- Browser-based Torsor
- WGSL shaders (same geometry, new backend)

**Multiphysics:**
- Heat transfer (analytical solutions for simple cases)
- Vibration modes (eigenvalue problems)
- Fluid flow (potential flow, analytical)

---

## Key Takeaways

### 1. Algebra Beats Iteration (When It Works)
For simple structural problems, analytical solutions (Roark's formulas) are:
- **Faster** (instant vs. hours of FEA)
- **Exact** (no mesh discretization error)
- **Parametric** (change dimension, instant recalculation)

FEA is still needed for complex geometries, but many engineering problems are **simpler than we think**.

### 2. Type Systems Are Your Friend
C++20 concepts let us write geometry code once and have it work with:
- Plain `double` (fast rendering)
- `Dual<double>` (automatic derivatives)
- Future `Gradient<double, N>` (optimization)

**The Pattern:**
```cpp
template <Scalar T>
T Sphere::eval(const std::array<T,3>& p) const {
    // Works for T=double, T=Dual<double>, T=Gradient<double,N>
}
```

### 3. Know When to Optimize (and When Not To)
**Don't optimize:** The renderer (finite differences are fine)
**Do optimize:** The solver (exact gradients matter)

**Engineering judgment** is knowing which epsilon errors you can tolerate.

### 4. Backend Separation Enables Evolution
By separating geometry (pure functions) from evaluation (CPU/GPU), we can:
- Develop on Mac (CPU, easy debugging)
- Deploy on Orin (GPU, fast execution)
- Add WebGPU later (browser, wide reach)

**Write once, run anywhere** isn't just for Java.

### 5. Documentation Is Part of the Product
We shipped:
- `README.md` (user-facing: what and how)
- `LESSONS.md` (developer-facing: why and future)

Six months from now, **you'll thank yourself** for writing this down.

---

## The Vision: MathCAD for 3D

**MathCAD** let you write `F = m * a` and see the answer update as you changed `m`.

**Torsor** lets you write `stress = M * c / I` and see the answer update as you change beam geometry.

**The Dream:**
```cpp
// Define a beam parametrically
ChannelBeam beam{
    .depth = design_var("depth", 6.0, 4.0, 12.0),  // Variable
    .length = 24.0
};

// Define constraints
constraint(max_stress < yield_strength);
constraint(deflection < span / 360);

// Solve
optimize(minimize(beam.weight()));
// → Torsor finds optimal depth = 6.47" (exact gradient descent)
```

**That's** the product. Everything we built tonight is the foundation.

---

## Phase 2: Assembly Design (v0.3 - Sunday Evening)

**Date:** Sunday, November 23, 2025
**Goal:** Extend Torsor from geometric modeling to component-level assembly design

### The Shift: From Geometry to Load Paths

**Saturday Night (v0.1-0.2):** We solved `f(x,y,z)` for geometry
**Sunday Evening (v0.3):** We solve the **Load Path**

```
Load on Shaft → Reaction on Bearing → Shear on Bolts → Bending on Channel
```

This is a **Sheaf of Physics**. The gluing condition is force equilibrium:
- Reaction force at bearing = Load applied to beam
- Bolt shear capacity > Bearing reaction
- Channel bending stress < Material yield

### Modular Refactoring

**Problem:** `main.cpp` was 726 lines of monolithic code
**Solution:** Extract into reusable headers

**Result:**
- `tensorcad_core.h` (107 lines) - Dual numbers, concepts
- `tensorcad_geometry.h` (295 lines) - Primitives, booleans, ChannelBeam
- `tensorcad_assembly.h` (62 lines) - Roark's formulas
- `tensorcad_render.h` (115 lines) - Raymarching, camera
- `main.cpp` (220 lines) - Application logic only

**Benefits:**
- **Reusability:** Both tools share geometry/math code
- **Maintainability:** Changes to physics don't affect rendering
- **Testability:** Each module can be verified independently
- **Evolution:** Easy to add WebGPU/Metal backends later

### The Assembly Designer

**Interactive CLI Tool** for real engineering workflows:

```bash
./assembly_designer

Enter radial load (lbs): 1500
Enter RPM: 1800
Enter torque (lb-in): 2500
Enter span (inches): 24
...

>>> PHASE 1: SHAFT & BEARING SELECTION
  [MATCH] Shaft Dia: 1.00" | Bearing: 7728T56 | Shaft SF: 2.3
  ✓ SELECTED: 7728T56 (Shaft Dia: 1.00")

>>> PHASE 2: C-CHANNEL SELECTION
  [MATCH] C 4 x 7.25 | Stress: 4250 psi | Weight: 7.25 lb/ft
  ✓ SELECTED: C 4 x 7.25

>>> PHASE 3: BOLTED CONNECTION
  Bolts: 2 x 1/2" Grade 5
  Shear Check: PASS
  Web Bearing Check: PASS

FINAL DESIGN SPECIFICATION:
1. SHAFT:      1.00" Dia (1045 Steel Q&T)
2. BEARINGS:   Sealmaster 7728T56
3. STRUCTURE:  C 4 x 7.25
4. FASTENERS:  2 x 1/2" Grade 5 Hex Bolts
```

### NASA TM-87354: Shaft Fatigue Analysis

**The Problem:** Shafts fail due to **cyclic loading** (fatigue), not static overload.

**Traditional Approach:** Conservative guesswork ("make it 2x bigger")

**Torsor Approach:** Exact fatigue life calculation

**Implementation:**
```cpp
ShaftResult analyze_shaft(double diameter, double torque, double moment,
                          double yield, double ultimate) {
    // 1. Marin Factors (NASA TM-87354)
    double k_surface = 0.85;       // Machined finish
    double k_size = (d < 2.0) ? 0.85 : 0.75;
    double k_reliability = 0.814;  // 99% reliability

    // 2. Endurance Limit
    double Se = 0.5 * ultimate * k_surface * k_size * k_reliability;

    // 3. Stresses
    double sigma_a = 32 * M / (π * d³);  // Alternating bending
    double tau_m = 16 * T / (π * d³);    // Mean torsion

    // 4. Modified Goodman Criterion
    double SF = 1 / (sigma_a/Se + √3*tau_m/Sut);

    return {SF, Se, sigma_a, SF >= 1.5};
}
```

**Why This Matters:**
- Accounts for **surface finish** (machined vs ground)
- Accounts for **size effect** (larger shafts are weaker per unit area)
- Accounts for **reliability** (99% vs 50% survival)
- Uses **Modified Goodman** for combined bending+torsion

**No guesswork. Pure physics.**

### Component Databases

Real engineering catalogs hardcoded for instant lookup:

**AISC C-Channels:**
```cpp
{"C 6 x 13", 6.0, 2.157, 0.343, 0.437, 3.83, 17.4, 5.80, 13.0}
{"C 5 x 9",  5.0, 1.885, 0.325, 0.320, 2.64, 8.90, 3.56, 9.0}
// depth, width, web_t, flange_t, area, Ix, Sx, weight
```

**Sealmaster Bearings (McMaster-Carr):**
```cpp
{"7728T56", 1.00, 2800, 6300, "Cast Iron"}
{"7728T71", 1.25, 5750, 5400, "Cast Iron"}
// part#, shaft_dia, load_cap, max_rpm, housing
```

**Future:** Load from CSV/JSON, query McMaster API, integrate vendor catalogs

### Design Iteration Workflow

**The User Experience:**

1. **Input Requirements** (load, RPM, torque, span)
2. **Tool Iterates** through component catalogs
3. **Physics Checks** at each level:
   - Bearing: Load capacity, RPM limit
   - Shaft: Fatigue SF >= 1.5 (NASA)
   - Channel: Bending stress < allowable (AISC)
   - Bolts: Shear capacity > reaction (SunCam)
4. **Output:** Optimized specification (lightest valid design)

**No Spreadsheets. No Manual Lookup. Instant.**

### Key Architectural Decisions

**1. Separation of Concerns**
- `assembly_designer.cpp` has NO geometry rendering
- `main.cpp` has NO component databases
- Shared physics in `tensorcad_assembly.h`

**2. Interactive vs Automated**
- Geometric tool: Hardcoded demo (visual output)
- Assembly tool: Interactive CLI (engineering spec output)
- Different use cases, different interfaces

**3. Makefile Build System**
```makefile
all: tensorcad assembly_designer

tensorcad: main.cpp $(HEADERS)
assembly_designer: assembly_designer.cpp

clean: rm -f tensorcad assembly_designer *.ppm
```

Both tools share headers, compile independently, run separately.

### Lessons from Phase 2

**1. Modular Architecture is Non-Negotiable**
Once you have 2 executables, monolithic code becomes unmaintainable. Extract early.

**2. Real Engineering Data Matters**
Toy examples are fun. Real McMaster part numbers and AISC tables make it **usable**.

**3. Interactive Beats Hardcoded**
Saturday's demo was impressive. Sunday's tool is **practical** - engineers can use it Monday.

**4. Physics Checkers are Composable**
```cpp
bool shaft_ok = analyze_shaft(...).passed;
bool bearing_ok = reaction < bearing.capacity;
bool channel_ok = stress < allowable;
bool bolt_ok = check_bolt_shear(...);

if (shaft_ok && bearing_ok && channel_ok && bolt_ok) {
    // Valid design!
}
```

Each check is independent. Easy to add more (deflection, vibration, buckling).

**5. The Load Path is a Sheaf**
Just like geometry (local sections + gluing = global shape), assemblies are:
- **Local sections:** Component specs (bearing capacity, channel properties)
- **Gluing maps:** Force equilibrium (reaction = load, bolt_shear > reaction)
- **Global section:** A valid assembly where all physics constraints hold

### Performance Notes

**Compilation Time:**
- Header-only templates mean longer compile time (0.5s → 1.2s)
- Acceptable for prototyping
- Future: Precompiled headers (PCH) for production

**Runtime:**
- Assembly designer: Instant (<0.1s for full iteration)
- No rendering, just arithmetic
- Could handle 1000+ component catalog in <1s

**Memory:**
- Component catalogs: ~1KB total (tiny)
- No dynamic allocation in hot loops
- Stack-based physics solvers

---

## Final Reflection

**What we proved (Phase 1 + 2):**
- Algebraic geometry (SDF) works for engineering CAD
- Analytical stress calculations eliminate FEA for simple cases
- Parametric design "just works" when geometry is pure functions
- C++20 concepts enable beautiful abstractions
- Modular architecture scales from demo to production
- Real engineering workflows need component databases
- You can build A LOT in one weekend

**What we learned:**
- Sheaf theory maps to both CSG operations AND load paths
- Dual numbers vs. Gradient types (renderer vs. solver)
- Register pressure matters (finite differences beat autodiff for rendering)
- Backend separation enables platform portability
- Interactive tools beat hardcoded demos for actual use
- Physics checkers are composable building blocks

**What's next (Phase 3 and beyond):**
- `Gradient<T,N>` for the parametric optimizer
- More component types (I-beams, motors, gearboxes)
- CSV/JSON catalog loading (McMaster API integration)
- GPU acceleration (Metal for Mac, CUDA for Orin)
- Web interface (WebGPU + WASM)
- The full vision: "Enter load, get optimized design" in real-time

---

**Author:** Built over one weekend (Saturday + Sunday) with Claude Code
**Timeline:**
- Saturday night: Geometric kernel, Roark's formulas, SDF rendering (v0.1-0.2)
- Sunday evening: Modular refactoring, assembly designer, NASA fatigue (v0.3)

**Inspiration:** MathCAD, sheaf theory, Roark's formulas, NASA TM-87354, industrial automation
**Philosophy:** *"Optimization can be replaced by algebra when the problem has the right symmetry."*

---

## Appendix: References

**Mathematical Foundation:**
- Roark's Formulas for Stress and Strain (9th edition)
- NASA TM-87354: Shaft Fatigue Life Under Combined Bending and Torsion
- AISC Steel Construction Manual (section properties)
- Sheaf Theory in Geometry and Logic (Mac Lane & Moerdijk)
- SunCam: Bolted Joint Design (bolt shear calculations)

**Technical Resources:**
- Inigo Quilez: SDF rendering techniques
- Signed Distance Field tutorials
- C++20 Concepts reference

**Related Work:**
- ImplicitCAD (Haskell-based SDF CAD)
- OpenSCAD (CSG-based parametric CAD)
- FreeCAD (B-Rep parametric CAD)

**What Makes Torsor Different:**
- Combines algebraic geometry with analytical stress analysis
- Uses sheaf theory as design principle
- Targets real-time optimization (future GPU solver)
- Clean C++20 implementation suitable for embedded systems (Orin)

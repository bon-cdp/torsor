# TensorCAD

**Parametric CAD + Analytical Stress Analysis**

A revolutionary approach to engineering CAD that combines algebraic geometry with analytical stress calculations, eliminating the need for iterative FEA for simple structural analysis.

## Concept

TensorCAD treats geometry as **algebraic functions** (Signed Distance Fields) rather than triangle meshes. This enables:

- **Parametric design**: Change a dimension, everything updates instantly
- **Analytical solutions**: Apply Roark's formulas directly without meshing
- **Differentiable geometry**: Automatic gradients via dual numbers
- **Joints as operations**: Welds are Unions, contacts are Intersections (Sheaf Theory)

## Features

### Geometry Primitives
- Sphere, Box, Cylinder, Beam (oriented rectangular prisms)
- C6x13 Channel Beam (AISC standard structural member)
- Boolean operations (Union, Intersection, Difference)

### Engineering Analysis
- **Section properties**: Area, moment of inertia, section modulus (calculated from algebraic definition)
- **Stress analysis**: Bending stress (σ = Mc/I), Hertzian contact stress
- **Safety factors**: Material yield strength vs. calculated stress
- **Roark's Formulas**: Chapter 8 (Beams), Chapter 14 (Contact)

### The "Sheaf Theory" Insight
Joints are **geometric intersections**, not separate entities:
- `Union(beam, support)` = Combined structure
- `Intersection(beam, support)` = Contact patch
- Stress distribution follows from SDF gradients

## Build & Run

```bash
# Compile
clang++ -std=c++20 -O3 main.cpp -o tensorcad

# Run
./tensorcad

# View output (ParaView or any PPM viewer)
open output.ppm
```

## Example Output

```
C6x13 American Standard Channel (24" long) on cylinder support
Point load: 1000 lbs
Max bending stress: 974.8 psi
Safety factor: 36.9 ✓ DESIGN OK
```

## Demo

The current demo shows:
1. **C6x13 channel beam** (24" long) resting on a **cylindrical support** (2" diameter)
2. **Section properties** calculated analytically from the beam geometry
3. **Stress analysis** using Roark's formulas for bending and contact
4. **3D render** of the combined structure

## Parametric Power

Change `beam.depth` from 6.0 to 8.0 in `main.cpp`:
- Section properties update **instantly**
- Stresses recalculate **analytically**
- No FEA mesh, no iteration, pure algebra

This is **MathCAD for 3D**.

## Tech Stack

- **Language**: C++20 (concepts, templates)
- **Math**: Dual numbers for automatic differentiation
- **Rendering**: Raymarching through SDFs
- **Output**: Binary PPM (ParaView compatible)

## Theory

Based on:
- **Implicit Geometry** (F-Rep / SDF)
- **Sheaf Theory** (local sections + gluing constraints)
- **Roark's Formulas for Stress and Strain** (9th edition)
- **AISC Steel Construction Manual** (section properties)

## Future

- More structural shapes (I-beams, angles, tees)
- More Roark's formulas (torsion, columns, shells)
- Optimization solver (minimize weight for target safety factor)
- STL import → algebraic fitting
- GPU acceleration (Metal/CUDA)

## Author

Built in one Saturday night by bon-cdp, combining:
- Industrial automation experience (Orin/Lexium/PyComm)
- Sheaf theory research (algebraic learning)
- Engineering knowledge (Roark's formulas)
- Modern C++ expertise (Google HEIR)

---

*"Optimization can be replaced by algebra when the problem has the right symmetry."*

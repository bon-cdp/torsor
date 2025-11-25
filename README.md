# Torsor

**Parametric CAD + Analytical Stress Analysis + Assembly Design**

A revolutionary approach to engineering CAD that combines algebraic geometry with analytical stress calculations, eliminating the need for iterative FEA for simple structural analysis.

## What's New in v0.3

🎉 **Modular Architecture**: Refactored into reusable header files (core, geometry, assembly, render)

🔧 **Assembly Design Tool**: Interactive CLI for selecting shaft/bearing/channel/bolt combinations

📐 **NASA TM-87354**: Shaft fatigue analysis with Marin factors and Modified Goodman criterion

🏗️ **AISC Integration**: Real component catalogs (C-channels, Sealmaster bearings)

## Concept

Torsor treats geometry as **algebraic functions** (Signed Distance Fields) rather than triangle meshes. This enables:

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

### Using Makefile (Recommended)

```bash
# Build both tools
make

# Run geometric modeling tool
./tensorcad

# Run assembly design tool
./assembly_designer

# View geometric output (ParaView or any PPM viewer)
open output.ppm
```

### Manual Compilation

```bash
# Geometric tool
clang++ -std=c++20 -O3 main.cpp -o tensorcad

# Assembly designer
clang++ -std=c++20 -O3 assembly_designer.cpp -o assembly_designer
```

## Tools

### 1. Geometric Modeling (`./tensorcad`)
- Parametric geometry using SDFs
- Analytical section properties
- Roark's stress formulas
- 3D raymarching visualization
- Outputs: PPM render

### 2. Assembly Designer (`./assembly_designer`)
- Interactive CLI for component selection
- Shaft-bearing-channel-bolt assemblies
- NASA TM-87354 shaft fatigue analysis
- AISC beam bending checks
- Bolt shear verification
- Outputs: Final specification report

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
- **Architecture**: Modular headers (core, geometry, assembly, render)
- **Math**: Dual numbers for automatic differentiation
- **Rendering**: Raymarching through SDFs
- **Output**: Binary PPM (ParaView compatible)
- **Build System**: Makefile with multiple executables

## File Structure

```
tensorcad/
├── tensorcad_core.h          # Dual numbers, concepts
├── tensorcad_geometry.h       # Primitives, boolean ops, ChannelBeam
├── tensorcad_assembly.h       # Roark's formulas
├── tensorcad_render.h         # Raymarching, camera
├── main.cpp                   # Geometric modeling tool
├── assembly_designer.cpp      # Assembly design tool
├── Makefile                   # Build system
└── README.md, LESSONS.md      # Documentation
```

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

Built by bon-cdp over one weekend:
- **v0.1-0.2** (Saturday night): Geometric kernel, Roark's formulas, SDF rendering
- **v0.3** (Sunday evening): Modular architecture, assembly design tool, NASA fatigue analysis

Combining:
- Industrial automation experience (Orin/Lexium/PyComm)
- Sheaf theory research (algebraic learning)
- Engineering knowledge (Roark's formulas, NASA TM-87354)
- Modern C++ expertise (Google HEIR)

---

*"Optimization can be replaced by algebra when the problem has the right symmetry."*

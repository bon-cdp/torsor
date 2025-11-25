// Torsor v0.3: Modular Parametric CAD with Engineering Analysis
// Geometric design tool using algebraic SDFs and Roark's formulas
// Compile: clang++ -std=c++20 -O3 main.cpp -o tensorcad

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// Torsor modular headers
#include "tensorcad_core.h"
#include "tensorcad_geometry.h"
#include "tensorcad_assembly.h"
#include "tensorcad_render.h"

// ============================================================================
// Output - PPM File Writer
// ============================================================================

void save_ppm(const std::string& filename, int width, int height,
              const std::vector<uint8_t>& buffer) {
    std::ofstream file(filename, std::ios::binary);

    // PPM P6 header (binary format for ParaView compatibility)
    file << "P6\n" << width << " " << height << "\n255\n";

    // Write raw RGB data (binary)
    file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());

    file.close();
    std::cout << "✓ Saved render to " << filename << " (ParaView compatible)" << std::endl;
}

// ============================================================================
// The Demo - Putting it all together
// ============================================================================

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Torsor v0.3 - Engineering Analysis Demo                ║\n";
    std::cout << "║  Parametric CAD + Roark's Stress Formulas                ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";

    // ========================================================================
    // PART A: Define the Geometry (Parametric!)
    // ========================================================================

    std::cout << "═══ GEOMETRY DEFINITION ═══\n\n";

    // C6x13 Channel Beam (24 inches long, horizontal)
    ChannelBeam beam{
        {0.0, 3.0, -10.0},        // Center position
        {1.0, 0.0, 0.0},          // Axis (along X)
        {0.0, 1.0, 0.0},          // Up (along Y)
        24.0                      // Length (inches)
    };

    std::cout << "✓ C6x13 American Standard Channel:\n";
    std::cout << "  - Length: " << beam.length << " inches\n";
    std::cout << "  - Depth: " << beam.depth << " inches\n";
    std::cout << "  - Flange width: " << beam.flange_width << " inches\n";
    std::cout << "  - Web thickness: " << beam.web_thickness << " inches\n\n";

    // Cylinder support (2" diameter, 12" long, vertical)
    Cylinder support{
        {0.0, 0.0, -10.0},        // Base center
        {0.0, 1.0, 0.0},          // Axis (vertical)
        1.0,                      // Radius (1 inch)
        6.0                       // Height (6 inches)
    };

    std::cout << "✓ Cylindrical Support:\n";
    std::cout << "  - Diameter: " << support.radius * 2 << " inches\n";
    std::cout << "  - Height: " << support.height << " inches\n\n";

    // Combined scene
    auto beam_geom = beam.get_geometry();
    Union<decltype(beam_geom), Cylinder> scene{beam_geom, support};

    // ========================================================================
    // PART B: Analytical Section Properties
    // ========================================================================

    std::cout << "═══ SECTION PROPERTIES (from algebraic definition) ═══\n\n";

    auto props = beam.calculate_properties();

    std::cout << "  Area (A):              " << props.area << " in²\n";
    std::cout << "  Moment of Inertia (Ix): " << props.Ix << " in⁴\n";
    std::cout << "  Moment of Inertia (Iy): " << props.Iy << " in⁴\n";
    std::cout << "  Section Modulus (Sx):   " << props.Sx << " in³\n";
    std::cout << "  Section Modulus (Sy):   " << props.Sy << " in³\n\n";

    // ========================================================================
    // PART C: Load Case & Stress Analysis (Roark's Formulas)
    // ========================================================================

    std::cout << "═══ STRESS ANALYSIS (Roark's Formulas) ═══\n\n";

    // Load case: Simply supported beam with center point load
    double point_load = 1000.0;  // 1000 lbs at center
    double span = beam.length;    // 24 inches

    std::cout << "Load Case: Simply supported, center point load\n";
    std::cout << "  Point load: " << point_load << " lbs\n";
    std::cout << "  Span: " << span << " inches\n\n";

    // Calculate maximum moment (Roark's Ch. 8)
    RoarksFormulas::SimpleBeamLoading loading{
        RoarksFormulas::SimpleBeamLoading::CENTER_POINT,
        point_load,
        span
    };
    double max_moment = loading.calculate_max_moment();

    std::cout << "  Max Bending Moment (M): " << max_moment << " lb-in\n";

    // Calculate bending stress (σ = Mc/I)
    RoarksFormulas::BeamStressAnalysis stress_calc{
        max_moment,
        beam.depth / 2.0,  // c = distance to extreme fiber
        props.Ix
    };
    double max_bending_stress = stress_calc.calculate_bending_stress();

    std::cout << "  Max Bending Stress (σ): " << max_bending_stress << " psi\n";

    // Contact stress at beam-cylinder interface (Roark's Ch. 14)
    RoarksFormulas::HertzianContact contact{
        point_load,              // Normal force
        30.0e6,                  // E for steel (30 Mpsi)
        1000.0,                  // Beam radius (approx as large)
        support.radius,          // Cylinder radius
        beam.web_thickness       // Contact length
    };
    double contact_stress = contact.calculate_max_contact_stress();

    std::cout << "  Contact Stress (Hertzian): " << contact_stress << " psi\n\n";

    // Safety factor (assuming structural steel, Fy = 36 ksi)
    double yield_strength = 36000.0;  // psi
    double safety_factor = yield_strength / max_bending_stress;

    std::cout << "  Material: A36 Structural Steel (Fy = 36 ksi)\n";
    std::cout << "  Safety Factor: " << safety_factor << "\n\n";

    if (safety_factor >= 1.5) {
        std::cout << "  ✓ DESIGN OK (SF >= 1.5)\n\n";
    } else {
        std::cout << "  ✗ DESIGN FAILS (SF < 1.5)\n\n";
    }

    // ========================================================================
    // PART D: The "Joint as Intersection" Concept
    // ========================================================================

    std::cout << "═══ SHEAF THEORY: JOINTS AS BOOLEAN OPERATIONS ═══\n\n";
    std::cout << "  Union (beam + support):      Combined structure\n";
    std::cout << "  Intersection (beam ∩ support): Contact patch\n";
    std::cout << "  → No separate 'weld' entity needed!\n";
    std::cout << "  → Contact area = geometric intersection\n";
    std::cout << "  → Stress distribution follows from SDF gradients\n\n";

    // ========================================================================
    // PART E: Render the Scene
    // ========================================================================

    std::cout << "═══ RENDERING ═══\n\n";

    Camera camera{
        {15.0, 5.0, 0.0},      // Position (to the side)
        {0.0, 3.0, -10.0},     // Look at beam center
        50.0                   // FOV
    };

    constexpr int width = 800;
    constexpr int height = 600;
    std::vector<uint8_t> image_buffer(width * height * 3);

    std::cout << "Rendering " << width << "x" << height << " scene...\n";

    for (int y = 0; y < height; y++) {
        if (y % 100 == 0) {
            std::cout << "  Progress: " << (y * 100 / height) << "%\n";
        }

        for (int x = 0; x < width; x++) {
            double screen_x = (2.0 * x / width) - 1.0;
            double screen_y = 1.0 - (2.0 * y / height);

            auto color = raymarch(scene, camera, screen_x, screen_y);

            size_t idx = (y * width + x) * 3;
            image_buffer[idx + 0] = color[0];
            image_buffer[idx + 1] = color[1];
            image_buffer[idx + 2] = color[2];
        }
    }

    std::cout << "  Progress: 100%\n\n";
    save_ppm("output.ppm", width, height, image_buffer);

    // ========================================================================
    // PART F: Parametric Demonstration
    // ========================================================================

    std::cout << "═══ PARAMETRIC POWER ═══\n\n";
    std::cout << "Try this experiment:\n";
    std::cout << "  1. Change beam.depth from 6.0 to 8.0 (line ~50)\n";
    std::cout << "  2. Recompile: clang++ -std=c++20 -O3 main.cpp -o tensorcad\n";
    std::cout << "  3. Run again: ./tensorcad\n";
    std::cout << "  → Section properties update INSTANTLY\n";
    std::cout << "  → Stresses recalculate ANALYTICALLY\n";
    std::cout << "  → No FEA mesh, no iteration, pure algebra!\n\n";

    std::cout << "This is MathCAD for 3D. This is Torsor.\n\n";

    return 0;
}

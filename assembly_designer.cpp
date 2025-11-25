// Torsor v0.3: Assembly Design Tool
// Interactive CLI for selecting shaft, bearing, channel, and bolt components
// Based on NASA TM-87354 (shaft fatigue) and AISC codes (beam bending)
// Compile: clang++ -std=c++20 -O3 assembly_designer.cpp -o assembly_designer

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>
#include <fstream>

// Torsor modular headers for visualization
#include "tensorcad_core.h"
#include "tensorcad_geometry.h"
#include "tensorcad_render.h"

const double PI = 3.14159265359;

// ============================================================================
// COMPONENT DATABASE
// ============================================================================

struct C_Channel {
    std::string name;
    double depth_in;
    double width_in;
    double web_thick_in;
    double flange_thick_in;
    double area_in2;
    double Ix_in4;
    double Sx_in3;  // Elastic Section Modulus
    double weight_lb_ft;
};

// AISC American Standard Channels (subset)
std::vector<C_Channel> get_channels() {
    return {
        {"C 6 x 13",   6.0, 2.157, 0.343, 0.437, 3.83, 17.4, 5.80, 13.0},
        {"C 6 x 10.5", 6.0, 2.034, 0.314, 0.326, 3.09, 15.2, 5.06, 10.5},
        {"C 5 x 9",    5.0, 1.885, 0.325, 0.320, 2.64, 8.90, 3.56, 9.0},
        {"C 4 x 7.25", 4.0, 1.721, 0.321, 0.296, 2.13, 4.59, 2.29, 7.25},
        {"C 3 x 6",    3.0, 1.596, 0.356, 0.273, 1.76, 2.07, 1.38, 6.0}
    };
}

struct MountedBearing {
    std::string part_number;
    double shaft_dia_in;
    double dynamic_load_cap_lb;
    double max_rpm;
    std::string housing_type;
};

// Sealmaster pillow block bearings (McMaster-Carr)
std::vector<MountedBearing> get_bearings() {
    return {
        {"7728T51", 0.50, 2600, 7300, "Cast Iron"},
        {"7728T52", 0.625, 2600, 7300, "Cast Iron"},
        {"7728T53", 0.75, 2600, 7300, "Cast Iron"},
        {"7728T56", 1.00, 2800, 6300, "Cast Iron"},
        {"7728T57", 1.125, 4350, 5400, "Cast Iron"},
        {"7728T71", 1.25, 5750, 5400, "Cast Iron"},
        {"7728T62", 1.50, 7300, 4100, "Cast Iron"},
        {"7728T65", 2.00, 9750, 3200, "Cast Iron"}
    };
}

struct Bolt {
    std::string grade;
    double diameter_in;
    double shear_strength_psi;    // Single shear
    double tensile_strength_psi;
};

// Common bolt grades
std::vector<Bolt> get_bolts() {
    return {
        {"Grade 5", 0.375, 50000, 120000},
        {"Grade 5", 0.50, 50000, 120000},
        {"Grade 5", 0.625, 50000, 120000},
        {"Grade 8", 0.375, 75000, 150000},
        {"Grade 8", 0.50, 75000, 150000},
        {"Grade 8", 0.625, 75000, 150000}
    };
}

// ============================================================================
// PHYSICS SOLVERS
// ============================================================================

// NASA TM-87354 Shaft Fatigue Analysis
struct ShaftResult {
    double safety_factor;
    double endurance_limit;
    double max_stress;
    bool passed;
};

ShaftResult analyze_shaft(double diameter, double torque_lb_in, double moment_lb_in,
                          double material_yield_psi, double material_ultimate_psi) {

    // 1. Marin Factors (surface finish, size, reliability)
    double k_surface = 0.85;       // Machined finish
    double k_size = (diameter < 2.0) ? 0.85 : 0.75;
    double k_reliability = 0.814;  // 99% reliability

    // 2. Endurance Limit (Se)
    double Se_prime = 0.5 * material_ultimate_psi;
    double Se = k_surface * k_size * k_reliability * Se_prime;

    // 3. Stresses (alternating bending + mean torsion)
    double sigma_a = (32.0 * moment_lb_in) / (PI * std::pow(diameter, 3));
    double tau_m = (16.0 * torque_lb_in) / (PI * std::pow(diameter, 3));

    // 4. Modified Goodman Failure Criterion
    double term1 = sigma_a / Se;
    double term2 = (std::sqrt(3.0) * tau_m) / material_ultimate_psi;

    double safety_factor = 1.0 / (term1 + term2);

    return {safety_factor, Se, sigma_a, safety_factor >= 1.5};
}

// Bolt Shear Check (SunCam method)
bool check_bolt_shear(double load_lb, int num_bolts, double bolt_dia, double shear_strength_psi) {
    double area = PI * std::pow(bolt_dia/2.0, 2);
    double total_capacity = num_bolts * area * shear_strength_psi;
    return total_capacity > (load_lb * 1.5);  // 1.5 SF required
}

// AISC Beam Bending Check
bool check_beam_bending(double moment_lb_in, double Sx_in3, double allowable_stress_psi) {
    double bending_stress = moment_lb_in / Sx_in3;
    return bending_stress < allowable_stress_psi;
}

// ============================================================================
// INTERACTIVE DESIGN TOOL
// ============================================================================

void print_header() {
    std::cout << "╔═══════════════════════════════════════════════════════════╗\n";
    std::cout << "║    Torsor v0.3 - ASSEMBLY DESIGN TOOL                    ║\n";
    std::cout << "║    Shaft-Bearing-Channel Selection & Analysis            ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════╝\n\n";
}

struct DesignInputs {
    double load_lb;
    double rpm;
    double torque_lb_in;
    double span_in;
    double load_location_in;
    double shaft_yield_psi;
    double shaft_ultimate_psi;
    double beam_allowable_stress_psi;
};

DesignInputs get_user_inputs() {
    DesignInputs inputs;

    std::cout << "═══ DESIGN REQUIREMENTS ═══\n\n";

    std::cout << "Enter radial load on shaft (lbs): ";
    std::cin >> inputs.load_lb;

    std::cout << "Enter operating speed (RPM): ";
    std::cin >> inputs.rpm;

    std::cout << "Enter torque (lb-in): ";
    std::cin >> inputs.torque_lb_in;

    std::cout << "Enter span between bearings (inches): ";
    std::cin >> inputs.span_in;

    std::cout << "Enter load location from left bearing (inches): ";
    std::cin >> inputs.load_location_in;

    std::cout << "\n═══ MATERIAL PROPERTIES ═══\n\n";

    std::cout << "Shaft Material (1045 Steel Q&T):\n";
    std::cout << "  Yield strength (psi) [default: 90000]: ";
    std::string input;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::getline(std::cin, input);
    inputs.shaft_yield_psi = input.empty() ? 90000.0 : std::stod(input);

    std::cout << "  Ultimate strength (psi) [default: 120000]: ";
    std::getline(std::cin, input);
    inputs.shaft_ultimate_psi = input.empty() ? 120000.0 : std::stod(input);

    std::cout << "\nBeam Material (A36 Steel):\n";
    std::cout << "  Allowable bending stress (psi) [default: 24000]: ";
    std::getline(std::cin, input);
    inputs.beam_allowable_stress_psi = input.empty() ? 24000.0 : std::stod(input);

    std::cout << "\n";
    return inputs;
}

void run_design_study(const DesignInputs& inputs) {
    std::cout << "════════════════════════════════════════════════════════════\n";
    std::cout << "DESIGN ANALYSIS\n";
    std::cout << "════════════════════════════════════════════════════════════\n\n";

    // --- PHASE 1: SHAFT & BEARING SELECTION ---
    std::cout << ">>> PHASE 1: SHAFT & BEARING SELECTION\n\n";

    std::vector<MountedBearing> valid_bearings;
    auto catalog_bearings = get_bearings();

    for (const auto& bearing : catalog_bearings) {
        // A. Check RPM limit
        if (inputs.rpm > bearing.max_rpm) continue;

        // B. Check bearing load capacity
        // Reaction forces (simply supported beam)
        double a = inputs.load_location_in;
        double b = inputs.span_in - a;
        double R1 = inputs.load_lb * b / inputs.span_in;
        double R2 = inputs.load_lb * a / inputs.span_in;
        double max_reaction = std::max(R1, R2);

        if (max_reaction > bearing.dynamic_load_cap_lb) continue;

        // C. Check shaft fatigue
        // Max moment at load point
        double max_moment = R1 * a;

        ShaftResult shaft_res = analyze_shaft(
            bearing.shaft_dia_in,
            inputs.torque_lb_in,
            max_moment,
            inputs.shaft_yield_psi,
            inputs.shaft_ultimate_psi
        );

        if (shaft_res.passed) {
            std::cout << "  [MATCH] Shaft Dia: " << std::setw(6) << std::setprecision(3) << bearing.shaft_dia_in
                      << "\" | Bearing: " << bearing.part_number
                      << " | Shaft SF: " << std::setprecision(2) << shaft_res.safety_factor << "\n";
            valid_bearings.push_back(bearing);
        }
    }

    if (valid_bearings.empty()) {
        std::cout << "\n  [FAIL] No valid shaft/bearing combinations found.\n";
        std::cout << "  Try reducing load, torque, or RPM.\n\n";
        return;
    }

    // Select smallest valid bearing
    MountedBearing selected_bearing = valid_bearings[0];
    std::cout << "\n  ✓ SELECTED: " << selected_bearing.part_number
              << " (Shaft Dia: " << selected_bearing.shaft_dia_in << "\")\n\n";


    // --- PHASE 2: STRUCTURAL SUPPORT SELECTION ---
    std::cout << ">>> PHASE 2: STRUCTURAL SUPPORT (C-CHANNEL) SELECTION\n\n";

    std::vector<C_Channel> valid_channels;
    auto catalog_channels = get_channels();

    // Reaction load on channel (bearing mounting)
    double a = inputs.load_location_in;
    double b = inputs.span_in - a;
    double bearing_load = inputs.load_lb * b / inputs.span_in;

    // Assume channel is simply supported over some span (e.g., 36")
    double channel_span = 36.0;
    double channel_moment = (bearing_load * channel_span) / 4.0;  // Center point load

    for (const auto& chan : catalog_channels) {
        double bending_stress = channel_moment / chan.Sx_in3;

        if (bending_stress < inputs.beam_allowable_stress_psi) {
            std::cout << "  [MATCH] " << std::setw(12) << std::left << chan.name
                      << " | Stress: " << std::setw(6) << std::right << static_cast<int>(bending_stress)
                      << " psi | Weight: " << chan.weight_lb_ft << " lb/ft\n";
            valid_channels.push_back(chan);
        }
    }

    if (valid_channels.empty()) {
        std::cout << "\n  [FAIL] No channels strong enough.\n\n";
        return;
    }

    // Select lightest valid channel
    C_Channel selected_channel = valid_channels.back();
    std::cout << "\n  ✓ SELECTED: " << selected_channel.name
              << " (" << selected_channel.weight_lb_ft << " lb/ft)\n\n";


    // --- PHASE 3: BOLTED CONNECTION CHECK ---
    std::cout << ">>> PHASE 3: BOLTED CONNECTION (BEARING TO CHANNEL)\n\n";

    int num_bolts = 2;
    double bolt_dia = 0.5;  // 1/2"
    double bolt_shear_strength = 50000.0;  // Grade 5

    // Side load (assume 20% of radial for conservatism)
    double side_load = bearing_load * 0.2;

    bool shear_ok = check_bolt_shear(side_load, num_bolts, bolt_dia, bolt_shear_strength);

    // Bearing stress on channel web
    double channel_bearing_stress = side_load / (num_bolts * bolt_dia * selected_channel.web_thick_in);
    double allow_bearing = 1.35 * 36000.0;  // AISC Fp

    std::cout << "  Bolts: " << num_bolts << " x 1/2\" Grade 5\n";
    std::cout << "  Shear Check: " << (shear_ok ? "PASS" : "FAIL") << "\n";
    std::cout << "  Web Bearing Stress: " << static_cast<int>(channel_bearing_stress)
              << " psi (Limit: " << static_cast<int>(allow_bearing) << " psi)\n";
    std::cout << "  Web Bearing Check: " << (channel_bearing_stress < allow_bearing ? "PASS" : "FAIL") << "\n\n";


    // --- FINAL SPECIFICATION ---
    std::cout << "════════════════════════════════════════════════════════════\n";
    std::cout << "FINAL DESIGN SPECIFICATION\n";
    std::cout << "════════════════════════════════════════════════════════════\n\n";

    std::cout << "1. SHAFT:      " << selected_bearing.shaft_dia_in << "\" Dia (1045 Steel Q&T)\n";
    std::cout << "2. BEARINGS:   Sealmaster " << selected_bearing.part_number
              << " (" << selected_bearing.housing_type << ")\n";
    std::cout << "3. STRUCTURE:  " << selected_channel.name
              << " (Web: " << selected_channel.web_thick_in << "\")\n";
    std::cout << "4. FASTENERS:  " << num_bolts << " x 1/2\" Grade 5 Hex Bolts\n\n";

    std::cout << "════════════════════════════════════════════════════════════\n\n";


    // --- VISUALIZATION ---
    std::cout << ">>> RENDERING ASSEMBLY FOR PARAVIEW\n\n";

    // Build the assembly geometry
    // 1. Shaft (horizontal cylinder)
    Cylinder shaft{
        {-inputs.span_in/2.0, 0.0, 0.0},   // Left end
        {1.0, 0.0, 0.0},                    // Horizontal axis
        selected_bearing.shaft_dia_in / 2.0, // Radius
        inputs.span_in                       // Length
    };

    // 2. Left bearing (slightly larger cylinder)
    double bearing_height = 2.0;  // Assume 2" tall bearing housing
    Cylinder bearing_left{
        {-inputs.span_in/2.0, 0.0, -bearing_height/2.0},
        {0.0, 0.0, 1.0},  // Vertical
        selected_bearing.shaft_dia_in * 1.2,  // Slightly larger than shaft
        bearing_height
    };

    // 3. Right bearing
    Cylinder bearing_right{
        {inputs.span_in/2.0, 0.0, -bearing_height/2.0},
        {0.0, 0.0, 1.0},  // Vertical
        selected_bearing.shaft_dia_in * 1.2,
        bearing_height
    };

    // 4. C-Channel support (under left bearing)
    ChannelBeam channel{
        {-inputs.span_in/2.0, 0.0, -bearing_height - selected_channel.depth_in/2.0},
        {1.0, 0.0, 0.0},  // Horizontal
        {0.0, 0.0, 1.0},  // Up
        channel_span       // 36" long
    };
    channel.depth = selected_channel.depth_in;
    channel.flange_width = selected_channel.width_in;
    channel.web_thickness = selected_channel.web_thick_in;
    channel.flange_thickness = selected_channel.flange_thick_in;

    // Combine all components
    auto bearing_pair = Union<Cylinder, Cylinder>{bearing_left, bearing_right};
    auto shaft_bearings = Union<Cylinder, decltype(bearing_pair)>{shaft, bearing_pair};
    auto channel_geom = channel.get_geometry();
    auto assembly = Union<decltype(shaft_bearings), decltype(channel_geom)>{shaft_bearings, channel_geom};

    // Multi-view rendering (engineering drawing style)
    constexpr int view_width = 400;
    constexpr int view_height = 400;
    constexpr int composite_width = view_width * 2;
    constexpr int composite_height = view_height * 2;

    std::vector<uint8_t> composite(composite_width * composite_height * 3, 255); // White background

    std::cout << "Rendering 4-view engineering drawing...\n";

    // Define 4 cameras for orthographic views
    std::vector<std::pair<std::string, Camera>> views = {
        {"TOP", make_ortho_camera("TOP", inputs.span_in)},
        {"ISO", make_ortho_camera("ISO", inputs.span_in)},
        {"FRONT", make_ortho_camera("FRONT", inputs.span_in)},
        {"RIGHT", make_ortho_camera("RIGHT", inputs.span_in)}
    };

    std::vector<std::pair<int, int>> positions = {
        {0, 0},                    // Top-left: TOP
        {view_width, 0},           // Top-right: ISO
        {0, view_height},          // Bottom-left: FRONT
        {view_width, view_height}  // Bottom-right: RIGHT
    };

    // Render each view
    for (size_t view_idx = 0; view_idx < views.size(); view_idx++) {
        auto& [view_name, cam] = views[view_idx];
        std::cout << "  Rendering " << view_name << " view...\n";

        std::vector<uint8_t> view_buffer(view_width * view_height * 3);

        for (int y = 0; y < view_height; y++) {
            for (int x = 0; x < view_width; x++) {
                double screen_x = (2.0 * x / view_width) - 1.0;
                double screen_y = 1.0 - (2.0 * y / view_height);

                auto color = raymarch(assembly, cam, screen_x, screen_y);

                size_t idx = (y * view_width + x) * 3;
                view_buffer[idx + 0] = color[0];
                view_buffer[idx + 1] = color[1];
                view_buffer[idx + 2] = color[2];
            }
        }

        // Copy view to composite at correct position
        auto [offset_x, offset_y] = positions[view_idx];
        for (int y = 0; y < view_height; y++) {
            for (int x = 0; x < view_width; x++) {
                int src_idx = (y * view_width + x) * 3;
                int dst_idx = ((offset_y + y) * composite_width + (offset_x + x)) * 3;
                composite[dst_idx + 0] = view_buffer[src_idx + 0];
                composite[dst_idx + 1] = view_buffer[src_idx + 1];
                composite[dst_idx + 2] = view_buffer[src_idx + 2];
            }
        }
    }

    // Save as PNG
    save_png("assembly_3view.png", composite_width, composite_height, composite);
    std::cout << "\n✓ Saved: assembly_3view.png\n\n";

    // Display in terminal with chafa
    std::cout << "═══ TERMINAL PREVIEW ═══\n\n";
    system("chafa --size 80x40 --colors full assembly_3view.png");

    // Auto-open in Preview.app
    std::cout << "\n✓ Opening in Preview.app...\n";
    std::cout << "  (Press Cmd+Tab to see high-quality view)\n\n";
    system("open -a Preview assembly_3view.png &");

    std::cout << "✓ Assembly design complete!\n";
    std::cout << "  Text spec: see above\n";
    std::cout << "  3D views:  assembly_3view.png (4-view engineering drawing)\n\n";
}

int main() {
    print_header();

    std::cout << "This tool automatically selects compatible components for a\n";
    std::cout << "shaft assembly (shaft + bearings + mounting channel + bolts)\n";
    std::cout << "based on your load requirements.\n\n";

    char cont = 'y';
    while (cont == 'y' || cont == 'Y') {
        DesignInputs inputs = get_user_inputs();
        run_design_study(inputs);

        std::cout << "Design another assembly? (y/n): ";
        std::cin >> cont;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "\n";
    }

    std::cout << "Thank you for using Torsor Assembly Designer!\n";
    std::cout << "For geometric modeling, use: ./tensorcad\n\n";

    return 0;
}

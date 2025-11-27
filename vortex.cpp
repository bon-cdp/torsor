// Torsor v0.1: Vortex Shedding Calculator
// Terminal-based interactive tool for analyzing vortex-induced vibration
// in HSS round members
//
// Compile: clang++ -std=c++20 -O3 vortex.cpp -o vortex \
//          -I/opt/homebrew/opt/ftxui/include \
//          -L/opt/homebrew/opt/ftxui/lib \
//          -lftxui-component -lftxui-dom -lftxui-screen

#include "vortex_physics.h"

#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/color.hpp>

#include <iostream>
#include <sstream>
#include <iomanip>

using namespace ftxui;

// ============================================================================
// COLOR MAPPING SYSTEM
// ============================================================================

enum class ColorVariable {
    UDL,      // Uniform distributed load (kN/m)
    MOMENT,   // Bending moment (kN·m)
    STRESS    // Bending stress (MPa)
};

// Map normalized value [0,1] to color (Blue → Green → Yellow → Red)
inline Color value_to_color(double normalized) {
    if (normalized < 0.2) return Color::Blue;
    else if (normalized < 0.4) return Color::Cyan;
    else if (normalized < 0.6) return Color::Green;
    else if (normalized < 0.8) return Color::Yellow;
    else return Color::Red;
}

// Get value at normalized position x along beam for given color variable
inline double get_value_at_position(const ModeResult& mode, double x_norm, ColorVariable var) {
    if (mode.distribution.x_m.empty()) return 0.0;

    size_t idx = static_cast<size_t>(x_norm * (mode.distribution.x_m.size() - 1));
    if (idx >= mode.distribution.x_m.size()) idx = mode.distribution.x_m.size() - 1;

    switch (var) {
        case ColorVariable::UDL:
            return mode.distribution.udl_kN_m;  // Constant along beam
        case ColorVariable::MOMENT:
            return mode.distribution.moment_kNm[idx];
        case ColorVariable::STRESS:
            return mode.distribution.stress_MPa[idx];
        default:
            return 0.0;
    }
}

// Get maximum value for color variable across all modes
inline double get_max_value(const std::vector<ModeResult>& modes, ColorVariable var) {
    double max_val = 0.0;
    for (const auto& mode : modes) {
        switch (var) {
            case ColorVariable::UDL:
                max_val = std::max(max_val, mode.distribution.udl_kN_m);
                break;
            case ColorVariable::MOMENT:
                max_val = std::max(max_val, mode.distribution.max_moment_kNm);
                break;
            case ColorVariable::STRESS:
                max_val = std::max(max_val, mode.distribution.max_stress_MPa);
                break;
        }
    }
    return max_val;
}

// ============================================================================
// APPLICATION STATE
// ============================================================================

struct AppState {
    // Input parameters (METRIC UNITS)
    int material_index = 0;
    double D_mm = 200.0;         // Outer diameter (mm)
    double t_mm = 6.0;           // Wall thickness (mm)
    double L_m = 6.0;            // Length (m)
    double damping = 0.01;       // 1% damping (typical for steel)
    double axial_force_kN = 0.0; // Axial force (kN, positive = tension)

    // Custom material (if selected)
    double custom_E_GPa = 200.0;       // Young's modulus (GPa)
    double custom_rho = 7850.0;        // Density (kg/m³)

    // UI state
    int active_input = 0;             // Which input field is active
    ColorVariable active_color_var = ColorVariable::STRESS;  // Start with stress

    // Results
    HSS_Member member;
    std::vector<VortexResults> results;

    // Update member and recalculate
    void update() {
        // Set material
        if (static_cast<size_t>(material_index) < MATERIALS.size() - 1) {
            member.material = MATERIALS[material_index];
        } else {
            member.material = {"Custom", custom_E_GPa, custom_rho};
        }

        // Set geometry (metric)
        member.D_mm = D_mm;
        member.t_mm = t_mm;
        member.L_m = L_m;
        member.damping_ratio = damping;

        // Calculate derived properties
        member.calculate_properties();

        // Run vortex analysis with axial force
        results = VortexPhysics::analyze_member(member, 2, axial_force_kN);
    }
};

// ============================================================================
// UI COMPONENTS
// ============================================================================

// Top panel: Logo and equations
Element render_header() {
    return vbox({
        text("╔═══════════════════════════════════════════════════════════════╗") | color(Color::Cyan),
        text("║                      Torsor  \\|/                              ║") | color(Color::Cyan) | bold,
        text("║              Vortex Shedding Calculator v0.1                  ║") | color(Color::Cyan),
        text("╚═══════════════════════════════════════════════════════════════╝") | color(Color::Cyan),
        text(""),
        hbox({
            text("ω = A√(EI/μL⁴)  ") | color(Color::Yellow),
            text("│") | color(Color::GrayDark),
            text("  V = f(ω,D,S)  ") | color(Color::Yellow),
            text("│") | color(Color::GrayDark),
            text("  F = C₁/(√(L/D)·√(ζ-C₂ρD²/M))·qₕ·D") | color(Color::Yellow)
        }) | center,
        separator()
    });
}

// Render a single mode shape curve using Braille characters
Element render_mode_shape(const std::string& bc_name, int mode, double max_deflection,
                         int width, int height, Color curve_color) {
    auto canvas = Canvas(width * 2, height * 4);  // Braille gives 2x4 resolution

    // Draw original member line (horizontal)
    int baseline_y = height * 4 / 2;
    for (int x = 0; x < width * 2; x++) {
        canvas.DrawPointLine(x, baseline_y, x, baseline_y, Color::GrayDark);
    }

    // Draw mode shape
    for (int x = 0; x < width * 2 - 1; x++) {
        double x_norm = static_cast<double>(x) / (width * 2);
        double x_next = static_cast<double>(x + 1) / (width * 2);

        double y1 = VortexPhysics::modal_shape(x_norm, 1.0, bc_name, mode);
        double y2 = VortexPhysics::modal_shape(x_next, 1.0, bc_name, mode);

        // Scale to canvas height (±height/4 from baseline)
        int canvas_y1 = baseline_y - static_cast<int>(y1 * height);
        int canvas_y2 = baseline_y - static_cast<int>(y2 * height);

        canvas.DrawPointLine(x, canvas_y1, x + 1, canvas_y2, curve_color);
    }

    return ftxui::canvas(std::move(canvas));
}

// Middle panel: Modal shapes with auto-scaling and color mapping
Element render_visualizations(const AppState& state) {
    if (state.results.empty()) {
        return text("Calculating...") | center;
    }

    std::vector<Element> bc_panels;

    // Get max value for color mapping across ALL boundary conditions
    double global_max_value = 0.0;
    for (const auto& result : state.results) {
        double bc_max = get_max_value(result.modes, state.active_color_var);
        global_max_value = std::max(global_max_value, bc_max);
    }

    for (const auto& result : state.results) {
        std::vector<Element> mode_shapes;

        // Title
        mode_shapes.push_back(text(result.bc_name) | bold | color(Color::Cyan) | center);

        // Calculate auto-scaling factor for this BC
        double max_deflection = 0.0;
        for (const auto& mode_res : result.modes) {
            for (int x = 0; x <= 100; x++) {
                double x_norm = x / 100.0;
                double deflection = VortexPhysics::modal_shape(x_norm, 1.0, result.bc_name, mode_res.mode_number);
                max_deflection = std::max(max_deflection, std::abs(deflection));
            }
        }

        // Canvas setup
        constexpr int canvas_width = 30;
        constexpr int canvas_height = 10;
        auto combined_canvas = Canvas(canvas_width * 2, canvas_height * 4);
        int baseline_y = canvas_height * 4 / 2;

        // Scale factor to fit modes in canvas (use 80% of height)
        double scale_factor = max_deflection > 0 ? (canvas_height * 4 * 0.4) / max_deflection : 1.0;

        // Draw baseline
        for (int x = 0; x < canvas_width * 2; x++) {
            combined_canvas.DrawPointLine(x, baseline_y, x, baseline_y, Color::GrayDark);
        }

        // Draw modes with color mapping
        for (const auto& mode_res : result.modes) {
            int mode = mode_res.mode_number;

            for (int x = 0; x < canvas_width * 2 - 1; x++) {
                double x_norm = static_cast<double>(x) / (canvas_width * 2);
                double x_next = static_cast<double>(x + 1) / (canvas_width * 2);

                // Get modal deflections
                double y1 = VortexPhysics::modal_shape(x_norm, 1.0, result.bc_name, mode);
                double y2 = VortexPhysics::modal_shape(x_next, 1.0, result.bc_name, mode);

                // Scale to canvas
                int canvas_y1 = baseline_y - static_cast<int>(y1 * scale_factor);
                int canvas_y2 = baseline_y - static_cast<int>(y2 * scale_factor);

                // Get color based on value at this position
                double value = get_value_at_position(mode_res, x_norm, state.active_color_var);
                double normalized = global_max_value > 0 ? value / global_max_value : 0.0;
                Color line_color = value_to_color(normalized);

                // Draw colored line
                combined_canvas.DrawPointLine(x, canvas_y1, x + 1, canvas_y2, line_color);
            }
        }

        mode_shapes.push_back(ftxui::canvas(std::move(combined_canvas)));

        // Mode info with max deflection and scale
        std::ostringstream scale_info;
        scale_info << "Max: " << std::fixed << std::setprecision(4) << max_deflection
                   << " | Scale: " << std::fixed << std::setprecision(1) << scale_factor;
        mode_shapes.push_back(text(scale_info.str()));  // Default color

        // Mode legends with intermediate values (split into multiple lines for readability)
        for (const auto& mode_res : result.modes) {
            // Line 1: Basic frequency and wind speed
            std::ostringstream oss1;
            oss1 << "Mode " << mode_res.mode_number << ": "
                 << "f=" << std::fixed << std::setprecision(2) << mode_res.freq_hz << " Hz"
                 << " | ω=" << std::fixed << std::setprecision(2) << mode_res.omega << " rad/s"
                 << " | V=" << std::fixed << std::setprecision(3) << mode_res.V_critical_ms << " m/s";
            mode_shapes.push_back(text(oss1.str()));  // Default color

            // Line 2: Force, stress, and Strouhal number
            std::ostringstream oss2;
            oss2 << "  F=" << std::fixed << std::setprecision(4) << mode_res.F_static_kN << " kN"
                 << " | σ_max=" << std::fixed << std::setprecision(2) << mode_res.distribution.max_stress_MPa << " MPa"
                 << " | S=" << std::fixed << std::setprecision(3) << mode_res.strouhal;
            mode_shapes.push_back(text(oss2.str()));  // Default color

            // Line 3: Amplification and damping parameters
            std::ostringstream oss3;
            oss3 << "  H=" << std::fixed << std::setprecision(1) << mode_res.amplification
                 << " | ωD²=" << std::fixed << std::setprecision(4) << mode_res.omega_D2 << " m²rad/s"
                 << " | ζ_eff=" << std::fixed << std::setprecision(5) << mode_res.damping_param;
            mode_shapes.push_back(text(oss3.str()));  // Default color
        }

        bc_panels.push_back(vbox(mode_shapes) | border);
    }

    // Arrange BC panels in 2 rows of 3 for better visibility
    std::vector<Element> top_row, bottom_row;
    for (size_t i = 0; i < bc_panels.size(); i++) {
        if (i < 3) {
            top_row.push_back(bc_panels[i]);
        } else {
            bottom_row.push_back(bc_panels[i]);
        }
    }

    return vbox({
        hbox(top_row) | center,
        hbox(bottom_row) | center
    });
}

// Color legend panel
Element render_color_legend(const AppState& state) {
    // Get variable name and units
    std::string var_name, units;
    double max_value = 0.0;

    // Get max across all results
    for (const auto& result : state.results) {
        double bc_max = get_max_value(result.modes, state.active_color_var);
        max_value = std::max(max_value, bc_max);
    }

    switch (state.active_color_var) {
        case ColorVariable::UDL:
            var_name = "UDL";
            units = "kN/m";
            break;
        case ColorVariable::MOMENT:
            var_name = "MOMENT";
            units = "kN·m";
            break;
        case ColorVariable::STRESS:
            var_name = "STRESS";
            units = "MPa";
            break;
    }

    // Create gradient bar
    std::vector<Element> gradient;
    for (int i = 0; i < 20; i++) {
        double norm = static_cast<double>(i) / 19.0;
        Color c = value_to_color(norm);
        gradient.push_back(text("█") | color(c));
    }

    // Format max value
    std::ostringstream max_str;
    max_str << std::fixed << std::setprecision(3) << max_value;

    return vbox({
        hbox({
            text("Showing: ") | color(Color::GrayLight),
            text(var_name) | bold | color(Color::Yellow),
            text(" (" + units + ")") | color(Color::GrayLight),
            text("  [Space to cycle]") | color(Color::GrayDark),
            text("  │  ") | color(Color::GrayDark),
            text("MAX = ") | color(Color::GrayLight),
            text(max_str.str() + " " + units) | bold | color(Color::Red)
        }) | center,
        hbox(gradient) | center,
        hbox({
            text("0.0") | color(Color::Blue),
            text(" ───────────── ") | color(Color::GrayDark),
            text(max_str.str()) | color(Color::Red)
        }) | center
    }) | border;
}

// Bottom panel: Inputs and outputs
Element render_inputs(const AppState& state) {
    auto material_name = (static_cast<size_t>(state.material_index) < MATERIALS.size())
                        ? MATERIALS[state.material_index].name
                        : "Custom";

    std::ostringstream oss_D, oss_t, oss_L, oss_damp, oss_axial;
    oss_D << std::fixed << std::setprecision(1) << state.D_mm;
    oss_t << std::fixed << std::setprecision(2) << state.t_mm;
    oss_L << std::fixed << std::setprecision(1) << state.L_m;
    oss_damp << std::fixed << std::setprecision(4) << state.damping;
    oss_axial << std::fixed << std::setprecision(1) << state.axial_force_kN;

    auto highlight_if_active = [&](int index, const std::string& label, const std::string& value) {
        auto elem = hbox({
            text(label + ": ") | color(Color::Yellow),
            text(value) | color(Color::White)
        });
        if (state.active_input == index) {
            return elem | bold | bgcolor(Color::Blue);
        }
        return elem;
    };

    return vbox({
        text("═══ INPUTS (Tab to cycle, ↑↓ to adjust) ═══") | color(Color::Cyan) | center,
        text(""),
        hbox({
            vbox({
                highlight_if_active(0, "Material", material_name),
                highlight_if_active(1, "Diameter (D)", oss_D.str() + " mm"),
                highlight_if_active(2, "Thickness (t)", oss_t.str() + " mm"),
            }) | flex,
            separator(),
            vbox({
                highlight_if_active(3, "Length (L)", oss_L.str() + " m"),
                highlight_if_active(4, "Damping (ζ)", oss_damp.str()),
                highlight_if_active(5, "Axial Force (F)", oss_axial.str() + " kN"),
            }) | flex,
            separator(),
            vbox({
                text("Calculated:") | color(Color::GrayLight),
                [&]() {
                    std::ostringstream oss;
                    oss << "  I = " << std::fixed << std::setprecision(0) << (state.member.I_m4 * 1e12) << " mm⁴";
                    return text(oss.str()) | color(Color::GrayLight);
                }(),
                [&]() {
                    std::ostringstream oss;
                    oss << "  A = " << std::fixed << std::setprecision(1) << (state.member.A_m2 * 1e6) << " mm²";
                    return text(oss.str()) | color(Color::GrayLight);
                }(),
                [&]() {
                    std::ostringstream oss;
                    oss << "  μ = " << std::fixed << std::setprecision(2) << state.member.mu_kg_m << " kg/m";
                    return text(oss.str()) | color(Color::GrayLight);
                }(),
            }) | flex
        }),
        separator(),
        text("q=quit  Space=color  Tab=input  ↑↓=adjust") | color(Color::GrayDark) | center
    }) | border;
}

// ============================================================================
// MAIN APPLICATION
// ============================================================================

int main() {
    auto screen = ScreenInteractive::Fullscreen();
    AppState state;

    // Initial calculation
    state.update();

    // Component for capturing keyboard input
    auto component = Renderer([&] {
        return vbox({
            render_header(),
            separator(),
            render_visualizations(state) | flex,
            separator(),
            render_color_legend(state),
            separator(),
            render_inputs(state)
        });
    });

    // Event handler
    component = CatchEvent(component, [&](Event event) {
        // Quit
        if (event == Event::Character('q') || event == Event::Character('Q')) {
            screen.ExitLoopClosure()();
            return true;
        }

        // Tab: cycle active input
        if (event == Event::Tab) {
            state.active_input = (state.active_input + 1) % 6;  // 6 inputs now
            return true;
        }

        // Space: cycle color variable (UDL → Moment → Stress)
        if (event == Event::Character(' ')) {
            switch (state.active_color_var) {
                case ColorVariable::UDL:
                    state.active_color_var = ColorVariable::MOMENT;
                    break;
                case ColorVariable::MOMENT:
                    state.active_color_var = ColorVariable::STRESS;
                    break;
                case ColorVariable::STRESS:
                    state.active_color_var = ColorVariable::UDL;
                    break;
            }
            return true;
        }

        // Arrow keys: adjust values
        bool changed = false;

        if (event == Event::ArrowUp) {
            switch (state.active_input) {
                case 0: state.material_index = (state.material_index + 1) % MATERIALS.size(); break;
                case 1: state.D_mm += 10.0; break;     // Increment by 10mm
                case 2: state.t_mm += 0.5; break;      // Increment by 0.5mm
                case 3: state.L_m += 0.5; break;       // Increment by 0.5m
                case 4: state.damping += 0.001; break;
                case 5: state.axial_force_kN += 10.0; break;  // Increment by 10kN
            }
            changed = true;
        }

        if (event == Event::ArrowDown) {
            switch (state.active_input) {
                case 0: state.material_index = (state.material_index - 1 + MATERIALS.size()) % MATERIALS.size(); break;
                case 1: if (state.D_mm > 20.0) state.D_mm -= 10.0; break;
                case 2: if (state.t_mm > 1.0) state.t_mm -= 0.5; break;
                case 3: if (state.L_m > 1.0) state.L_m -= 0.5; break;
                case 4: if (state.damping > 0.001) state.damping -= 0.001; break;
                case 5: state.axial_force_kN -= 10.0; break;  // Can be negative (compression)
            }
            changed = true;
        }

        // Recalculate if inputs changed
        if (changed) {
            state.update();
            return true;
        }

        return false;
    });

    screen.Loop(component);

    std::cout << "\nThank you for using Torsor!\n";
    std::cout << "Keyboard-driven engineering tools for the modern era.\n\n";

    return 0;
}

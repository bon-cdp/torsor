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
// APPLICATION STATE
// ============================================================================

struct AppState {
    // Input parameters
    int material_index = 0;
    double D = 8.0;              // inches
    double t = 0.25;             // inches
    double L = 240.0;            // inches (20 feet)
    double damping = 0.01;       // 1% damping (typical for steel)

    // Custom material (if selected)
    double custom_E = 29.0e6;
    double custom_rho = 0.284;

    // UI state
    int active_input = 0;        // Which input field is active
    int color_variable = 0;      // 0=deflection, 1=moment, 2=stress, 3=UDL

    // Results
    HSS_Member member;
    std::vector<VortexResults> results;

    // Update member and recalculate
    void update() {
        // Set material
        if (static_cast<size_t>(material_index) < MATERIALS.size() - 1) {
            member.material = MATERIALS[material_index];
        } else {
            member.material = {"Custom", custom_E, custom_rho};
        }

        // Set geometry
        member.D = D;
        member.t = t;
        member.L = L;
        member.damping_ratio = damping;

        // Calculate derived properties
        member.calculate_properties();

        // Run vortex analysis
        results = VortexPhysics::analyze_member(member, 2);
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

// Middle panel: Modal shapes for all boundary conditions
Element render_visualizations(const AppState& state) {
    if (state.results.empty()) {
        return text("Calculating...") | center;
    }

    // Color for each mode
    std::vector<Color> mode_colors = {Color::Green, Color::Blue, Color::Magenta};

    std::vector<Element> bc_panels;

    for (const auto& result : state.results) {
        std::vector<Element> mode_shapes;

        // Title
        mode_shapes.push_back(text(result.bc_name) | bold | color(Color::Cyan) | center);

        // Draw each mode overlaid
        auto combined_canvas = Canvas(30 * 2, 10 * 4);
        int baseline_y = 10 * 4 / 2;

        // Draw baseline
        for (int x = 0; x < 30 * 2; x++) {
            combined_canvas.DrawPointLine(x, baseline_y, x, baseline_y, Color::GrayDark);
        }

        // Overlay modes
        for (size_t mode_idx = 0; mode_idx < result.modes.size(); mode_idx++) {
            int mode = result.modes[mode_idx].mode_number;
            Color mode_color = mode_colors[mode_idx % mode_colors.size()];

            for (int x = 0; x < 30 * 2 - 1; x++) {
                double x_norm = static_cast<double>(x) / (30 * 2);
                double x_next = static_cast<double>(x + 1) / (30 * 2);

                double y1 = VortexPhysics::modal_shape(x_norm, 1.0, result.bc_name, mode);
                double y2 = VortexPhysics::modal_shape(x_next, 1.0, result.bc_name, mode);

                int canvas_y1 = baseline_y - static_cast<int>(y1 * 8);
                int canvas_y2 = baseline_y - static_cast<int>(y2 * 8);

                combined_canvas.DrawPointLine(x, canvas_y1, x + 1, canvas_y2, mode_color);
            }
        }

        mode_shapes.push_back(ftxui::canvas(std::move(combined_canvas)));

        // Mode legends
        for (size_t i = 0; i < result.modes.size(); i++) {
            const auto& mode_res = result.modes[i];
            Color mode_color = mode_colors[i % mode_colors.size()];

            std::ostringstream oss;
            oss << "Mode " << mode_res.mode_number << ": "
                << std::fixed << std::setprecision(1) << mode_res.freq_hz << " Hz, "
                << std::fixed << std::setprecision(0) << (mode_res.V_critical * 0.0568182) << " mph";

            mode_shapes.push_back(text(oss.str()) | color(mode_color));
        }

        bc_panels.push_back(vbox(mode_shapes) | border);
    }

    return hbox(bc_panels) | center;
}

// Bottom panel: Inputs and outputs
Element render_inputs(const AppState& state) {
    auto material_name = (static_cast<size_t>(state.material_index) < MATERIALS.size())
                        ? MATERIALS[state.material_index].name
                        : "Custom";

    std::ostringstream oss_D, oss_t, oss_L, oss_damp;
    oss_D << std::fixed << std::setprecision(3) << state.D;
    oss_t << std::fixed << std::setprecision(3) << state.t;
    oss_L << std::fixed << std::setprecision(1) << state.L;
    oss_damp << std::fixed << std::setprecision(4) << state.damping;

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
        text("═══ INPUTS (Tab to cycle, ↑↓ to adjust, Shift to change color mapping) ═══") | color(Color::Cyan) | center,
        text(""),
        hbox({
            vbox({
                highlight_if_active(0, "Material", material_name),
                highlight_if_active(1, "Diameter (D)", oss_D.str() + " in"),
                highlight_if_active(2, "Thickness (t)", oss_t.str() + " in"),
            }) | flex,
            separator(),
            vbox({
                highlight_if_active(3, "Length (L)", oss_L.str() + " in"),
                highlight_if_active(4, "Damping (ζ)", oss_damp.str()),
                text(""),
            }) | flex,
            separator(),
            vbox({
                text("Calculated:") | color(Color::GrayLight),
                text("  I = " + std::to_string(state.member.I).substr(0, 6) + " in⁴") | color(Color::GrayLight),
                text("  A = " + std::to_string(state.member.A).substr(0, 6) + " in²") | color(Color::GrayLight),
                text("  μ = " + std::to_string(state.member.mu).substr(0, 6) + " lb/in") | color(Color::GrayLight),
            }) | flex
        }),
        separator(),
        text("q=quit  h=help  s=save") | color(Color::GrayDark) | center
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
            state.active_input = (state.active_input + 1) % 5;
            return true;
        }

        // Shift: cycle color variable
        if (event == Event::Character(' ')) {  // Using space as shift substitute for now
            state.color_variable = (state.color_variable + 1) % 4;
            return true;
        }

        // Arrow keys: adjust values
        bool changed = false;

        if (event == Event::ArrowUp) {
            switch (state.active_input) {
                case 0: state.material_index = (state.material_index + 1) % MATERIALS.size(); break;
                case 1: state.D += 0.5; break;
                case 2: state.t += 0.01; break;
                case 3: state.L += 12.0; break;
                case 4: state.damping += 0.001; break;
            }
            changed = true;
        }

        if (event == Event::ArrowDown) {
            switch (state.active_input) {
                case 0: state.material_index = (state.material_index - 1 + MATERIALS.size()) % MATERIALS.size(); break;
                case 1: if (state.D > 1.0) state.D -= 0.5; break;
                case 2: if (state.t > 0.01) state.t -= 0.01; break;
                case 3: if (state.L > 12.0) state.L -= 12.0; break;
                case 4: if (state.damping > 0.001) state.damping -= 0.001; break;
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

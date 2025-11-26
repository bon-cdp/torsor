#ifndef VORTEX_PHYSICS_H
#define VORTEX_PHYSICS_H

// Torsor v0.1: Vortex Shedding Analysis - Physics Engine
// Calculates natural frequencies, critical wind speeds, and dynamic forces
// for HSS round members under vortex-induced vibration

#include <cmath>
#include <string>
#include <vector>
#include <array>

const double PI = 3.14159265358979323846;

// ============================================================================
// MATERIAL DATABASE
// ============================================================================

struct Material {
    std::string name;
    double E_psi;              // Young's modulus (psi)
    double density_lb_in3;     // Density (lb/in³)
};

// Standard structural steel materials
const std::vector<Material> MATERIALS = {
    {"A36 Steel", 29.0e6, 0.284},
    {"A572-50 Steel", 29.0e6, 0.284},
    {"A500 Grade B HSS", 29.0e6, 0.284},
    {"Custom", 0.0, 0.0}  // User-defined
};

// ============================================================================
// BOUNDARY CONDITION COEFFICIENTS
// ============================================================================

struct BoundaryCondition {
    std::string name;
    std::vector<double> A_coeffs;  // A coefficients for each mode (1-indexed)

    // TODO: Complete table from reference
    // Currently using provided values for modes 1 and 2
};

const std::vector<BoundaryCondition> BOUNDARY_CONDITIONS = {
    {"Fixed-Hinged", {0.0, 15.4, 50.0}},      // Mode 0 (unused), Mode 1, Mode 2
    {"Fixed-Fixed",  {0.0, 22.4, 61.7}},
    {"Hinged-Hinged", {0.0, 9.87, 39.5}}
};

// ============================================================================
// HSS MEMBER PROPERTIES
// ============================================================================

struct HSS_Member {
    // Geometric properties
    double D;              // Outer diameter (in)
    double t;              // Wall thickness (in)
    double L;              // Length (in)

    // Material properties
    Material material;

    // Derived properties (calculated)
    double I;              // Moment of inertia (in⁴)
    double A;              // Cross-sectional area (in²)
    double mu;             // Mass per unit length (lb·s²/in² = lbm/in in consistent units)
    double r_gyration;     // Radius of gyration (in)

    // Dynamic properties
    double damping_ratio;  // Damping ratio (ζ, dimensionless, typically 0.001-0.02 for steel)

    // Calculate derived properties
    void calculate_properties() {
        // Area: A = π(D²-(D-2t)²)/4 = π(D² - D² + 4Dt - 4t²)/4 = πt(D-t)
        double D_inner = D - 2.0 * t;
        A = PI * (D*D - D_inner*D_inner) / 4.0;  // in²

        // Moment of inertia: I = π(D⁴-(D-2t)⁴)/64
        I = PI * (std::pow(D, 4) - std::pow(D_inner, 4)) / 64.0;  // in⁴

        // Mass per unit length: μ = ρ * A (lb/in)
        // Note: In FPS system with g = 386.4 in/s², to get consistent units for dynamics:
        // μ in dynamics equation has units lb·s²/in (divide by g to get from lbf to lbm)
        mu = material.density_lb_in3 * A;  // lb/in (mass per unit length)

        // Radius of gyration: r = sqrt(I/A)
        r_gyration = std::sqrt(I / A);  // in
    }
};

// ============================================================================
// VORTEX SHEDDING RESULTS
// ============================================================================

struct ModeResult {
    int mode_number;
    double omega;              // Natural frequency (rad/s)
    double freq_hz;            // Natural frequency (Hz)
    double V_critical;         // Critical wind speed (in/s, convert to mph for display)
    double F_static;           // Static force approximation (lb)
    double omega_axial;        // Natural frequency with axial load (rad/s)
    double amplification;      // Dynamic amplification factor |H(f)|

    // Intermediate calculations (for visualization)
    double omega_D2;           // ωD² parameter for Strouhal number selection
    double qh;                 // Velocity pressure (psi)
    double alpha;              // Axial load parameter
};

struct VortexResults {
    std::string bc_name;       // Boundary condition name
    std::vector<ModeResult> modes;  // Results for each mode
};

// ============================================================================
// PHYSICS CALCULATIONS
// ============================================================================

namespace VortexPhysics {

// Calculate natural frequency for a given mode
// ω = A * sqrt(E * I / (μ * L⁴))
// Units: [rad/s] = [1] * sqrt([psi]*[in⁴] / ([lb/in]*[in⁴]))
//                = sqrt([lb/in²]*[in⁴] / ([lb/in]*[in⁴]))
//                = sqrt([lb*in³] / [lb*in⁵])
//                = sqrt(1/[in²]) = 1/[in] × sqrt([in²/s²]) → need to check!
// Actually: E in psi, μ needs to be in (lb·s²/in)/in = lb·s²/in²
inline double calculate_natural_frequency(double E, double I, double mu, double L, double A_coeff) {
    // Convert mu from lb/in to lb·s²/in² (divide by g = 386.4 in/s²)
    double mu_dynamic = mu / 386.4;  // lb·s²/in²

    // ω = A * sqrt(E*I / (mu*L⁴))
    double omega = A_coeff * std::sqrt(E * I / (mu_dynamic * std::pow(L, 4)));
    return omega;  // rad/s
}

// Calculate critical wind speed using Strouhal number approach
// Different formulas based on ωD² value
inline double calculate_critical_wind_speed(double omega, double D) {
    double omega_D2 = omega * D * D;  // (rad/s) * in² = in²·rad/s

    double V;  // in/s

    if (omega_D2 <= 0.5) {
        V = 6.0 * omega * D;
    } else if (omega_D2 < 0.75) {
        V = 3.0 * omega * D + 1.5 / D;
    } else {
        V = 5.0 * omega * D;
    }

    return V;  // in/s
}

// Calculate static force approximation from vortex shedding
// F = C1 / (sqrt(L/D) * (ζ - C2*(ρ*D²)/M)^0.5) * qh * D
// where qh = 0.6 * V²
inline double calculate_static_force(double V, double D, double L, double damping, double mu, double rho_air = 7.346e-5) {
    // Constants
    const double C1 = 3.0;
    const double C2 = 0.6;

    // Velocity pressure: qh = 0.6 * V² (assuming air density incorporated)
    // Standard: qh = 0.5 * ρ * V², but formula uses 0.6*V² empirically
    double qh = 0.6 * V * V;  // (in/s)² - needs density factor

    // Air density at sea level: approximately 0.0765 lb/ft³ = 4.43e-5 lb/in³
    // Actually using dynamic pressure in psi: qh = 0.5 * rho * V² where rho ~ 7.346e-5 lb·s²/in⁴
    qh = 0.5 * rho_air * V * V;  // psi

    // Denominator calculation
    double sqrt_L_D = std::sqrt(L / D);
    double mass_term = C2 * rho_air * D * D / mu;
    double denom_inner = damping - mass_term;

    // Handle case where damping is very low
    if (denom_inner <= 0) {
        denom_inner = 1e-6;  // Prevent division by zero
    }

    double denom = sqrt_L_D * std::sqrt(denom_inner);

    // Static force
    double F = (C1 / denom) * qh * D;  // lb

    return F;
}

// Calculate natural frequency with axial load effect
// ω_a = ω * sqrt(1 ± α²/n²)
// α = F*L²/(E*I*π²)
inline double calculate_axial_frequency(double omega, double F, double L, double E, double I, int mode) {
    double alpha_squared = (F * L * L) / (E * I * PI * PI);
    double mode_squared = mode * mode;

    // Use + for tension, - for compression
    // For vortex shedding, typically assume tension from wind drag
    double factor = 1.0 + alpha_squared / mode_squared;

    if (factor < 0) factor = 0;  // Can't have imaginary frequency

    double omega_axial = omega * std::sqrt(factor);
    return omega_axial;
}

// Calculate dynamic amplification factor for SDOF system
// |H(f)| = 1 / sqrt((1-r²)² + (2ζr)²)
// where r = f/fn (frequency ratio), typically use r = 1 for maximum
inline double calculate_amplification(double damping, double freq_ratio = 1.0) {
    double r = freq_ratio;
    double zeta = damping;

    double term1 = (1.0 - r*r) * (1.0 - r*r);
    double term2 = (2.0 * zeta * r) * (2.0 * zeta * r);

    double H = 1.0 / std::sqrt(term1 + term2);

    return H;
}

// TODO: Temporary value - will be replaced with proper calculation
const double TEMP_AMPLIFICATION = 50.0;

// Analyze a single mode
inline ModeResult analyze_mode(const HSS_Member& member, const BoundaryCondition& bc, int mode) {
    ModeResult result;
    result.mode_number = mode;

    // 1. Natural frequency
    result.omega = calculate_natural_frequency(
        member.material.E_psi,
        member.I,
        member.mu,
        member.L,
        bc.A_coeffs[mode]
    );
    result.freq_hz = result.omega / (2.0 * PI);

    // 2. Critical wind speed
    result.omega_D2 = result.omega * member.D * member.D;
    result.V_critical = calculate_critical_wind_speed(result.omega, member.D);

    // 3. Static force
    result.F_static = calculate_static_force(
        result.V_critical,
        member.D,
        member.L,
        member.damping_ratio,
        member.mu
    );

    // 4. Axial effect on frequency
    result.alpha = (result.F_static * member.L * member.L) / (member.material.E_psi * member.I * PI * PI);
    result.omega_axial = calculate_axial_frequency(
        result.omega,
        result.F_static,
        member.L,
        member.material.E_psi,
        member.I,
        mode
    );

    // 5. Dynamic amplification
    // TODO: Use proper calculation instead of temporary value
    result.amplification = TEMP_AMPLIFICATION;  // calculate_amplification(member.damping_ratio, 1.0);

    // Velocity pressure for reference
    result.qh = 0.5 * 7.346e-5 * result.V_critical * result.V_critical;

    return result;
}

// Analyze all boundary conditions and modes
inline std::vector<VortexResults> analyze_member(const HSS_Member& member, int num_modes = 2) {
    std::vector<VortexResults> all_results;

    for (const auto& bc : BOUNDARY_CONDITIONS) {
        VortexResults bc_results;
        bc_results.bc_name = bc.name;

        for (int mode = 1; mode <= num_modes && static_cast<size_t>(mode) < bc.A_coeffs.size(); mode++) {
            ModeResult mode_result = analyze_mode(member, bc, mode);
            bc_results.modes.push_back(mode_result);
        }

        all_results.push_back(bc_results);
    }

    return all_results;
}

// Modal shape function for visualization
// Returns deflection at position x along length L
// Mode shapes for different boundary conditions:
// Fixed-Fixed: y(x) = (1 - cos(2πnx/L))
// Fixed-Hinged: y(x) = sin(πnx/L) - sinh(βnx/L) (approximate)
// Hinged-Hinged: y(x) = sin(πnx/L)
inline double modal_shape(double x, double L, const std::string& bc_name, int mode) {
    double xi = x / L;  // Normalized position [0, 1]

    if (bc_name == "Hinged-Hinged") {
        return std::sin(mode * PI * xi);
    } else if (bc_name == "Fixed-Fixed") {
        // Approximate as (1 - cos(2πnξ))
        return 1.0 - std::cos(2.0 * PI * mode * xi);
    } else if (bc_name == "Fixed-Hinged") {
        // Simplified: use sine with slight modification
        return std::sin(mode * PI * xi) * (1.0 + 0.3 * xi);
    }

    return 0.0;
}

} // namespace VortexPhysics

#endif // VORTEX_PHYSICS_H

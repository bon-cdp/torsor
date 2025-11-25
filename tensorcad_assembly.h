#ifndef TENSORCAD_ASSEMBLY_H
#define TENSORCAD_ASSEMBLY_H

// Torsor Assembly: Engineering Analysis and Roark's Formulas
// Part of Torsor v0.3 - Modular Architecture

#include <cmath>

// ============================================================================
// Engineering Analysis - Roark's Formulas for Stress and Strain
// ============================================================================

namespace RoarksFormulas {
    // Chapter 8: Beams - Bending Stress
    // σ = M * c / I
    // where: M = bending moment, c = distance from neutral axis, I = moment of inertia
    struct BeamStressAnalysis {
        double moment;                  // Applied bending moment (lb-in)
        double max_fiber_distance;      // Distance from neutral axis to extreme fiber (in)
        double moment_of_inertia;       // Moment of inertia (in⁴)

        double calculate_bending_stress() const {
            return moment * max_fiber_distance / moment_of_inertia;  // psi
        }
    };

    // Chapter 8: Beams - Maximum Bending Moment for Simple Support
    // Simply supported beam with uniform load: M_max = w*L²/8
    // Simply supported beam with center point load: M_max = P*L/4
    struct SimpleBeamLoading {
        enum LoadType { UNIFORM, CENTER_POINT };
        LoadType type;
        double load;    // w (lb/in) for uniform, P (lb) for point
        double span;    // L (in)

        double calculate_max_moment() const {
            if (type == UNIFORM) {
                return load * span * span / 8.0;  // lb-in
            } else {
                return load * span / 4.0;  // lb-in
            }
        }
    };

    // Chapter 14: Bodies in Contact - Hertzian Contact Stress
    // Cylinder on cylinder (crossed cylinders)
    // σ_max = 0.591 * sqrt(P * E / (r1 * r2 * L))
    // where P = normal force, E = elastic modulus, r = radii, L = contact length
    struct HertzianContact {
        double force;                   // Normal force (lb)
        double elastic_modulus;         // E (psi), for steel ~30e6 psi
        double radius1, radius2;        // Radii of contacting cylinders (in)
        double contact_length;          // Length of contact (in)

        double calculate_max_contact_stress() const {
            // Simplified for steel-on-steel contact
            double E_eff = elastic_modulus;  // For same material
            return 0.591 * std::sqrt(force * E_eff / (radius1 * radius2 * contact_length));
        }
    };
}

#endif // TENSORCAD_ASSEMBLY_H

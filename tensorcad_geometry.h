#ifndef TENSORCAD_GEOMETRY_H
#define TENSORCAD_GEOMETRY_H

// Torsor Geometry: Primitives, Boolean Operations, and Engineering Shapes
// Part of Torsor v0.3 - Modular Architecture

#include "tensorcad_core.h"
#include <array>
#include <cmath>

// ============================================================================
// PART 1: Geometry Primitives - The "Local Sections" (Sheaf Theory)
// ============================================================================

struct Sphere {
    std::array<double, 3> center;
    double radius;

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // SDF of sphere: distance from center minus radius
        T dx = p[0] - center[0];
        T dy = p[1] - center[1];
        T dz = p[2] - center[2];
        T dist = sqrt(dx*dx + dy*dy + dz*dz);
        return dist - radius;
    }
};

struct Box {
    std::array<double, 3> corner;  // Min corner
    std::array<double, 3> size;    // Dimensions

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // SDF of axis-aligned box
        T dx = abs(p[0] - (corner[0] + size[0]/2)) - size[0]/2;
        T dy = abs(p[1] - (corner[1] + size[1]/2)) - size[1]/2;
        T dz = abs(p[2] - (corner[2] + size[2]/2)) - size[2]/2;

        // Distance to box surface
        T exterior = sqrt(max(dx, T(0))*max(dx, T(0)) +
                         max(dy, T(0))*max(dy, T(0)) +
                         max(dz, T(0))*max(dz, T(0)));
        T interior = min(max(dx, max(dy, dz)), T(0));
        return exterior + interior;
    }
};

struct Cylinder {
    std::array<double, 3> base;    // Bottom center
    std::array<double, 3> axis;    // Direction (should be normalized)
    double radius;
    double height;

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // Vector from base to point
        T dx = p[0] - base[0];
        T dy = p[1] - base[1];
        T dz = p[2] - base[2];

        // Project onto axis to get height along cylinder
        T h = dx*axis[0] + dy*axis[1] + dz*axis[2];

        // Perpendicular distance (radial component)
        T px = dx - h*axis[0];
        T py = dy - h*axis[1];
        T pz = dz - h*axis[2];
        T radial_dist = sqrt(px*px + py*py + pz*pz);

        // Distance from radius
        T dr = radial_dist - radius;

        // Distance from height bounds
        T dh = abs(h - height/2.0) - height/2.0;

        // Combine (inside cylinder when both negative)
        return max(dr, dh);
    }
};

struct Beam {
    std::array<double, 3> center;     // Center point
    std::array<double, 3> axis;       // Length direction (normalized)
    std::array<double, 3> up;         // Up direction (normalized, perpendicular to axis)
    double length;                    // Length along axis
    double width;                     // Width perpendicular to axis
    double height;                    // Height in the up direction

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // Transform point to beam's local coordinate system
        T dx = p[0] - center[0];
        T dy = p[1] - center[1];
        T dz = p[2] - center[2];

        // Project onto beam's axes
        T l = dx*axis[0] + dy*axis[1] + dz*axis[2];        // Length coordinate
        T h = dx*up[0] + dy*up[1] + dz*up[2];              // Height coordinate

        // Width coordinate (perpendicular to both axis and up)
        // Using cross product: width_dir = axis × up
        double wx = axis[1]*up[2] - axis[2]*up[1];
        double wy = axis[2]*up[0] - axis[0]*up[2];
        double wz = axis[0]*up[1] - axis[1]*up[0];
        T w = dx*wx + dy*wy + dz*wz;                       // Width coordinate

        // Distance from each dimension's bounds
        T dl = abs(l) - length/2.0;
        T dw = abs(w) - width/2.0;
        T dh = abs(h) - height/2.0;

        // Box SDF: max of all three distances for interior,
        // plus exterior distance for points outside
        T exterior = sqrt(max(dl, T(0))*max(dl, T(0)) +
                         max(dw, T(0))*max(dw, T(0)) +
                         max(dh, T(0))*max(dh, T(0)));
        T interior = min(max(dl, max(dw, dh)), T(0));
        return exterior + interior;
    }
};

// ============================================================================
// PART 2: Boolean Operations - The "Gluing Maps" (Sheaf Theory)
// ============================================================================

template <ImplicitField<double> A, ImplicitField<double> B>
struct Union {
    A shape_a;
    B shape_b;

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // Union is the minimum distance (closest surface wins)
        return min(shape_a.eval(p), shape_b.eval(p));
    }
};

template <ImplicitField<double> A, ImplicitField<double> B>
struct Intersection {
    A shape_a;
    B shape_b;

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // Intersection is the maximum distance (both must be inside)
        return max(shape_a.eval(p), shape_b.eval(p));
    }
};

template <ImplicitField<double> A, ImplicitField<double> B>
struct Difference {
    A shape_a;
    B shape_b;

    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        // Difference is A minus B
        return max(shape_a.eval(p), -shape_b.eval(p));
    }
};

// ============================================================================
// PART 3: Engineering Shapes - Standard Structural Members
// ============================================================================

// C6x13 American Standard Channel Beam
// Parametric definition allows instant recalculation when dimensions change
struct ChannelBeam {
    std::array<double, 3> center;     // Center point of the channel
    std::array<double, 3> axis;       // Length direction (normalized)
    std::array<double, 3> up;         // Up direction (web direction, normalized)
    double length;                    // Beam length

    // AISC C6x13 standard dimensions (inches)
    double depth = 6.0;               // Overall depth
    double flange_width = 2.157;      // Flange width
    double web_thickness = 0.343;     // Web thickness
    double flange_thickness = 0.437;  // Flange thickness

    // Construct the channel from 3 beams: web + 2 flanges
    Union<Beam, Union<Beam, Beam>> get_geometry() const {
        // Web (vertical center piece)
        Beam web{
            center,
            axis,
            up,
            length,
            web_thickness,
            depth
        };

        // Top flange (horizontal)
        double flange_offset = depth / 2.0 - flange_thickness / 2.0;
        std::array<double, 3> top_flange_center = {
            center[0] + up[0] * flange_offset,
            center[1] + up[1] * flange_offset,
            center[2] + up[2] * flange_offset
        };

        Beam top_flange{
            top_flange_center,
            axis,
            up,
            length,
            flange_width,
            flange_thickness
        };

        // Bottom flange (horizontal)
        std::array<double, 3> bot_flange_center = {
            center[0] - up[0] * flange_offset,
            center[1] - up[1] * flange_offset,
            center[2] - up[2] * flange_offset
        };

        Beam bot_flange{
            bot_flange_center,
            axis,
            up,
            length,
            flange_width,
            flange_thickness
        };

        // Combine all three beams
        Union<Beam, Beam> flanges{top_flange, bot_flange};
        return Union<Beam, Union<Beam, Beam>>{web, flanges};
    }

    // Analytical section properties (Roark's Chapter 8 / AISC Steel Manual)
    struct SectionProperties {
        double area;              // Cross-sectional area (in²)
        double Ix;                // Moment of inertia about x-axis (in⁴)
        double Iy;                // Moment of inertia about y-axis (in⁴)
        double Sx;                // Section modulus about x-axis (in³)
        double Sy;                // Section modulus about y-axis (in³)
    };

    SectionProperties calculate_properties() const {
        // For C6x13 (from AISC tables, but we calculate from geometry)
        double A = flange_width * flange_thickness * 2 +
                   (depth - 2*flange_thickness) * web_thickness;

        // Moment of inertia about strong axis (Ix)
        // Using parallel axis theorem for flanges
        double flange_y = (depth - flange_thickness) / 2.0;
        double Ix_flanges = 2 * (flange_width * flange_thickness * flange_thickness * flange_thickness / 12.0 +
                                 flange_width * flange_thickness * flange_y * flange_y);

        double web_height = depth - 2*flange_thickness;
        double Ix_web = web_thickness * web_height * web_height * web_height / 12.0;

        double Ix = Ix_flanges + Ix_web;

        // Moment of inertia about weak axis (Iy) - approximate
        double Iy = 2 * (flange_thickness * flange_width * flange_width * flange_width / 12.0) +
                    (web_height * web_thickness * web_thickness * web_thickness / 12.0);

        // Section modulus
        double Sx = Ix / (depth / 2.0);
        double Sy = Iy / (flange_width / 2.0);

        return {A, Ix, Iy, Sx, Sy};
    }

    // Template eval for SDF compatibility
    template <Scalar T>
    T eval(const std::array<T, 3>& p) const {
        auto geom = get_geometry();
        return geom.eval(p);
    }
};

#endif // TENSORCAD_GEOMETRY_H

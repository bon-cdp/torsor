#ifndef TENSORCAD_RENDER_H
#define TENSORCAD_RENDER_H

// Torsor Render: Raymarching and Visualization
// Part of Torsor v0.3 - Modular Architecture

#include "tensorcad_core.h"
#include <array>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <string>

// PNG support via stb_image_write
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// ============================================================================
// The Renderer - Raymarching through the Implicit Field
// ============================================================================

struct Camera {
    std::array<double, 3> position;
    std::array<double, 3> look_at;
    double fov;  // Field of view in degrees
};

// PNG output function
inline void save_png(const std::string& filename, int width, int height,
                     const std::vector<uint8_t>& buffer) {
    stbi_write_png(filename.c_str(), width, height, 3,
                   buffer.data(), width * 3);
}

// Create orthographic camera for engineering views
inline Camera make_ortho_camera(const std::string& view_name,
                                double assembly_size,
                                const std::array<double, 3>& center = {0,0,0}) {
    // Distance from object (far enough to avoid perspective distortion)
    double distance = assembly_size * 3.0;

    // Use small FOV to approximate orthographic projection
    double fov = 10.0;  // Narrow FOV = less perspective distortion

    Camera cam;
    cam.look_at = center;
    cam.fov = fov;

    if (view_name == "FRONT") {
        // Looking along +Y axis (from front)
        cam.position = {center[0], center[1] - distance, center[2]};

    } else if (view_name == "TOP") {
        // Looking down along -Z axis
        cam.position = {center[0], center[1], center[2] + distance};

    } else if (view_name == "RIGHT") {
        // Looking along -X axis (from right side)
        cam.position = {center[0] + distance, center[1], center[2]};

    } else if (view_name == "ISO") {
        // Isometric view (for context)
        cam.position = {
            center[0] + distance * 0.7,
            center[1] - distance * 0.7,
            center[2] + distance * 0.5
        };
        cam.fov = 45.0;  // Wider FOV for 3D view
    }

    return cam;
}

// Compute surface normal using the gradient (this is where Dual numbers shine!)
template <ImplicitField<double> G>
std::array<double, 3> compute_normal(const G& geometry, const std::array<double, 3>& p) {
    // We evaluate the geometry with Dual numbers to get gradients
    constexpr double epsilon = 0.0001;  // Smaller epsilon for more accurate normals

    // Evaluate with perturbations to get gradient via finite differences
    auto fx = geometry.eval(std::array{p[0] + epsilon, p[1], p[2]}) -
              geometry.eval(std::array{p[0] - epsilon, p[1], p[2]});
    auto fy = geometry.eval(std::array{p[0], p[1] + epsilon, p[2]}) -
              geometry.eval(std::array{p[0], p[1] - epsilon, p[2]});
    auto fz = geometry.eval(std::array{p[0], p[1], p[2] + epsilon}) -
              geometry.eval(std::array{p[0], p[1], p[2] - epsilon});

    // Normalize
    double len = std::sqrt(fx*fx + fy*fy + fz*fz);
    return {fx/len, fy/len, fz/len};
}

template <ImplicitField<double> G>
std::array<uint8_t, 3> raymarch(const G& geometry, const Camera& cam,
                                double screen_x, double screen_y) {
    // Build camera coordinate system
    double aspect = 800.0 / 600.0;
    double fov_rad = cam.fov * M_PI / 180.0;
    double h = std::tan(fov_rad / 2.0);

    // Forward vector: camera looking toward look_at
    std::array<double, 3> forward = {
        cam.look_at[0] - cam.position[0],
        cam.look_at[1] - cam.position[1],
        cam.look_at[2] - cam.position[2]
    };
    double forward_len = std::sqrt(forward[0]*forward[0] + forward[1]*forward[1] + forward[2]*forward[2]);
    forward = {forward[0]/forward_len, forward[1]/forward_len, forward[2]/forward_len};

    // World up vector (assume +Y is up)
    std::array<double, 3> world_up = {0.0, 1.0, 0.0};

    // Right vector: forward × world_up
    std::array<double, 3> right = {
        forward[1]*world_up[2] - forward[2]*world_up[1],
        forward[2]*world_up[0] - forward[0]*world_up[2],
        forward[0]*world_up[1] - forward[1]*world_up[0]
    };
    double right_len = std::sqrt(right[0]*right[0] + right[1]*right[1] + right[2]*right[2]);
    right = {right[0]/right_len, right[1]/right_len, right[2]/right_len};

    // Up vector: right × forward
    std::array<double, 3> up = {
        right[1]*forward[2] - right[2]*forward[1],
        right[2]*forward[0] - right[0]*forward[2],
        right[0]*forward[1] - right[1]*forward[0]
    };

    // Ray direction: forward + screen_x*right*aspect*h + screen_y*up*h
    std::array<double, 3> dir = {
        forward[0] + screen_x * right[0] * aspect * h + screen_y * up[0] * h,
        forward[1] + screen_x * right[1] * aspect * h + screen_y * up[1] * h,
        forward[2] + screen_x * right[2] * aspect * h + screen_y * up[2] * h
    };

    // Normalize direction
    double len = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir = {dir[0]/len, dir[1]/len, dir[2]/len};

    // Raymarch: start at camera, step along ray until we hit surface
    std::array<double, 3> pos = cam.position;
    constexpr double max_dist = 100.0;
    constexpr double threshold = 0.001;  // Smaller threshold for better hit detection
    double traveled = 0.0;

    for (int i = 0; i < 100; i++) {
        double dist = geometry.eval(pos);

        if (dist < threshold) {
            // Hit! Compute lighting based on normal
            auto normal = compute_normal(geometry, pos);

            // Simple diffuse lighting (light from top-right)
            std::array<double, 3> light_dir = {0.5, 0.5, -0.5};
            double light_len = std::sqrt(light_dir[0]*light_dir[0] +
                                        light_dir[1]*light_dir[1] +
                                        light_dir[2]*light_dir[2]);
            light_dir = {light_dir[0]/light_len, light_dir[1]/light_len, light_dir[2]/light_len};

            double diffuse = std::max(0.0,
                normal[0]*light_dir[0] +
                normal[1]*light_dir[1] +
                normal[2]*light_dir[2]);

            // Convert to color (blue-ish theme)
            uint8_t intensity = static_cast<uint8_t>(diffuse * 200 + 55);
            return {
                static_cast<uint8_t>(intensity * 0.5),
                static_cast<uint8_t>(intensity * 0.7),
                intensity
            };
        }

        // Step forward by the distance (safe because SDF gives minimum distance)
        traveled += dist;
        if (traveled > max_dist) break;

        pos = {
            pos[0] + dir[0] * dist,
            pos[1] + dir[1] * dist,
            pos[2] + dir[2] * dist
        };
    }

    // Miss: return background color (dark gray)
    return {20, 20, 30};
}

#endif // TENSORCAD_RENDER_H

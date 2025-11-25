#ifndef TENSORCAD_CORE_H
#define TENSORCAD_CORE_H

// Torsor Core: Dual Numbers and Type Concepts
// Part of Torsor v0.3 - Modular Architecture

#include <array>
#include <cmath>
#include <concepts>

// ============================================================================
// PART 1: The "Tensor" Foundation - Dual Numbers for Automatic Differentiation
// ============================================================================
// A Dual number carries both a value and its derivative.
// This implements forward-mode autodiff via operator overloading.

template <std::floating_point T>
struct Dual {
    T val;   // The actual value
    T grad;  // The derivative

    // Constructor for constants (derivative is zero)
    constexpr Dual(T v) : val(v), grad(T(0)) {}
    constexpr Dual(T v, T g) : val(v), grad(g) {}

    // Arithmetic operators implementing calculus rules
    constexpr Dual operator+(const Dual& rhs) const {
        return {val + rhs.val, grad + rhs.grad};
    }

    constexpr Dual operator-(const Dual& rhs) const {
        return {val - rhs.val, grad - rhs.grad};
    }

    constexpr Dual operator*(const Dual& rhs) const {
        // Product rule: (fg)' = f'g + fg'
        return {val * rhs.val, grad * rhs.val + val * rhs.grad};
    }

    constexpr Dual operator/(const Dual& rhs) const {
        // Quotient rule: (f/g)' = (f'g - fg') / g²
        return {val / rhs.val, (grad * rhs.val - val * rhs.grad) / (rhs.val * rhs.val)};
    }

    constexpr Dual operator-() const {
        return {-val, -grad};
    }

    // Comparison (based on value only, gradient doesn't affect ordering)
    constexpr bool operator<(const Dual& rhs) const { return val < rhs.val; }
    constexpr bool operator>(const Dual& rhs) const { return val > rhs.val; }
};

// Standard math functions extended to Dual numbers
template <std::floating_point T>
Dual<T> sqrt(const Dual<T>& x) {
    T sq = std::sqrt(x.val);
    return {sq, x.grad / (T(2) * sq)};  // d/dx sqrt(x) = 1/(2*sqrt(x))
}

template <std::floating_point T>
Dual<T> abs(const Dual<T>& x) {
    return {std::abs(x.val), x.val >= T(0) ? x.grad : -x.grad};
}

template <std::floating_point T>
Dual<T> min(const Dual<T>& a, const Dual<T>& b) {
    return (a.val < b.val) ? a : b;
}

template <std::floating_point T>
Dual<T> max(const Dual<T>& a, const Dual<T>& b) {
    return (a.val > b.val) ? a : b;
}

// Overloads for plain floating point (when not using Dual numbers)
template <std::floating_point T>
constexpr T min(T a, T b) {
    return (a < b) ? a : b;
}

template <std::floating_point T>
constexpr T max(T a, T b) {
    return (a > b) ? a : b;
}

// ============================================================================
// PART 2: The Concept - What makes a valid "Shape"?
// ============================================================================

// A Scalar is either a raw float/double OR a Dual number
template<typename T>
concept Scalar = std::floating_point<T> || requires(T t) {
    { t.val } -> std::floating_point;
    { t.grad } -> std::floating_point;
};

// An ImplicitField is anything that can be evaluated at a 3D point
// and returns a signed distance (negative = inside, positive = outside)
template<typename G, typename T>
concept ImplicitField = Scalar<T> && requires(const G& geometry, std::array<T, 3> point) {
    { geometry.eval(point) } -> std::convertible_to<T>;
};

#endif // TENSORCAD_CORE_H

#pragma once

#include <AMReX_REAL.H>

/**
 * @file CRESCENT_Constants.H
 * @brief Centralized constants for the spectral solver
 *
 * This file defines physical constants, numerical parameters, and filter
 * settings used throughout the CRESCENT solver stack. Using constexpr ensures these
 * are compile-time constants with zero runtime overhead.
 */

namespace crescent::constants {

// ============================================================================
// Dealiasing Filter Constants
// ============================================================================

/// 2/3 dealiasing cutoff ratio (Orszag's 2/3 rule for pseudospectral methods)
constexpr amrex::Real TWO_THIRDS_CUTOFF = 2.0 / 3.0;

/// 1/3 dealiasing cutoff (complementary)
constexpr amrex::Real ONE_THIRD_CUTOFF = 1.0 / 3.0;

// ============================================================================
// Hou-Li Exponential Filter Constants
// ============================================================================

/// Mild Hou-Li filter strength (for general use, L0, surface derivatives)
/// sigma = exp(-alpha * kappa^p)
constexpr amrex::Real HOULIGHER_FILTER_ALPHA_MILD = 33000.0;

/// Aggressive Hou-Li filter strength (for L1 nested AMR, DNO operator)
/// Used to suppress zigzag oscillations in refined regions
constexpr amrex::Real HOULIGHER_FILTER_ALPHA_AGGRESSIVE = 33000.0;

/// Hou-Li filter power (focuses filtering on high frequencies)
/// Higher p means sharper transition (more selective high-freq removal)
constexpr amrex::Real HOULIGHER_FILTER_POWER = 20.0;

// ============================================================================
// Mathematical Constants
// ============================================================================

/// Pi constant
constexpr amrex::Real PI = 3.14159265358979323846;

/// 2*Pi (for wavenumber calculations)
constexpr amrex::Real PI2 = 2.0 * PI;

// ============================================================================
// Physical Constants
// ============================================================================

/// Gravitational acceleration (m/s²)
constexpr amrex::Real GRAVITY = 9.81;

/// Inverse gravity (for efficiency in kernels)
constexpr amrex::Real INV_GRAVITY = 1.0 / GRAVITY;

// ============================================================================
// Numerical Precision
// ============================================================================

/// Small epsilon for floating-point comparisons
constexpr amrex::Real EPSILON = 1.0e-12;

/// Zero threshold for numerical operations
constexpr amrex::Real SMALL_NUMBER = 1.0e-14;

// ============================================================================
// FFT Normalization
// ============================================================================

/// AMReX forwardThenBackward gives N² scaling factor
/// This function computes the normalization factor for a given grid size
constexpr amrex::Real fftNormalization(int nx, int ny) {
    return 1.0 / static_cast<amrex::Real>(nx * ny);
}

} // namespace crescent::constants

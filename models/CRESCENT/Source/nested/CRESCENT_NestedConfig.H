#ifndef CRESCENT_NESTED_CONFIG_H
#define CRESCENT_NESTED_CONFIG_H

#include <AMReX_REAL.H>

namespace crescent::nested {

struct AlgorithmConfig
{
    int use_crescent_algorithm = 0;                   // 0=single-level L0 only, 1=enable CRESCENT for nested L1
    amrex::Real l1_box_x_lo = 0.0;                  // L1 refined box lower x boundary from input
    amrex::Real l1_box_x_hi = 0.0;                  // L1 refined box upper x boundary from input
    amrex::Real l1_box_y_lo = 0.0;                  // L1 refined box lower y boundary from input
    amrex::Real l1_box_y_hi = 0.0;                  // L1 refined box upper y boundary from input

    amrex::Real roi_accurate_fraction = 0.25;       // Margin for Ω_A (default: 25%)
    amrex::Real roi_spectral_fraction = 0.15;       // Margin for Ω_S (default: 15%)
    amrex::Real damping_steepness = 2.0;            // Damping-path shape parameter; ignored when C3 repair is active
    amrex::Real blending_steepness = 2.0;           // Accepted input; no active σ builder currently uses this knob
    amrex::Real blending_tensor_alpha = 1.0;        // Tensor-product σ profile parameter: 1=standard tensor form, larger values delay blend turn-on near Ω_A
    amrex::Real blending_superellipse_gamma = 1.0;  // Superellipse σ side-ramp shaping exponent: >1 delays blend turn-on near Ω_A
    amrex::Real blending_superellipse_p = 4.0;      // Superellipse σ exponent: larger values make the 2D contours more rectangle-like
    int blending_correction_y_smoothing_passes = 0; // Optional 1-2-1 passes on the blend correction field along y
    int state_interpolation_order = 2;              // L0->L1 eta/phi interpolation: 1=bilinear, 2=quadratic (default)
    int derivative_interpolation_order = 2;         // L0->L1 derivative interpolation: 1=bilinear, 2=quadratic (default)
    int blending_function_mode = 0;                 // 0=2D Euclidean, 1=sequential 1D, 2=harmonic 2D, 3=tensor-product 2D, 4=harmonic 2D + septic smoothing, 5=superellipse 2D
    int use_sequential_blending = 0;                // 0=2D Euclidean (default), 1=sequential 1D (preserves y-uniformity)

    int use_c3_repair = 1;                          // 0=use damping path, 1=use C3-repair (default: 1)
    int c3_apply_x = 1;                             // C3-repair direction switch: 0=skip x-pass, 1=apply x-pass (default)
    int c3_apply_y = 1;                             // C3-repair direction switch: 0=skip y-pass, 1=apply y-pass (default)
    int c3_use_interface = 1;                       // 0=cell-center, 1=interface-based (geometrically correct, default: 1)
    int c3_buffer_width = 24;                       // C3-repair buffer width in grid points (user override, default: auto-compute from Ω_A)
    amrex::Real c3_boundary_penalty = 0.001;        // C3-repair boundary magnitude constraint weight (default: 0.001, tuned)
    amrex::Real c3_coefficient_regularization = 1.0; // C3-repair coefficient magnitude control (L2 penalty on alpha, default: 1.0)
    int c3_x_pass_y_smoothing_passes = 0;           // Optional 1-2-1 smoothing passes on x-pass correction along y (default: 0)
    int c3_x_pass_corner_taper = 0;                 // Optional taper of x-pass inside top/bottom y-buffers to reduce corner coupling
    int c3_x_alpha_y_smoothing_passes = 0;          // Optional 1-2-1 smoothing passes on rowwise x-pass coefficients along y
    amrex::Real c3_x_alpha_y_regularization = 0.0;  // Coupled x-pass coefficient smoothness penalty along y (default: 0)
    amrex::Real c3_y_alpha_x_regularization = 0.0;  // Coupled y-pass coefficient smoothness penalty along x (default: 0)
    amrex::Real c3_x_alpha_y_mean_regularization = 0.0; // Coupled x-pass penalty on deviation from transverse-mean coefficient field
    amrex::Real c3_y_alpha_x_mean_regularization = 0.0; // Coupled y-pass penalty on deviation from transverse-mean coefficient field
    int c3_repair_after_filter = 0;                 // Optional second C3 repair after L1 Hou-Li filtering
    int c3_x_transverse_harmonics = 0;              // Periodic Fourier harmonics for x-pass coefficient field along y (0=linewise solve)
    int c3_y_transverse_harmonics = 0;              // Periodic Fourier harmonics for y-pass coefficient field along x (0=linewise solve)
    int c3_x_alpha_y_projection_harmonics = 0;      // Periodic Fourier projection of linewise x-pass coefficients along y
    int c3_y_alpha_x_projection_harmonics = 0;      // Periodic Fourier projection of linewise y-pass coefficients along x

    int filter_periodic_copies_before_linear_derivatives = 1; // 0=use post-C3 copies for linear eta/phi derivatives, 1=use filtered copies (default)
    int linear_derivative_spectral_filter_mode = 0;           // 0=2/3 cutoff only (default), 1=2/3+Hou-Li, 2=2/3+Hou-Li+low-pass
    int enable_post_timestep_low_pass_filter = 0;             // 0=preserve current post-timestep behavior (default), 1=apply low-pass filter after surf_update
    int enable_post_rk4_state_filter_and_blend = 1;           // 0=skip post-RK4 state overwrite, 1=run nested post-RK4 filter/blend (default)
    int enable_final_output_blended_derivative_filter = 0;    // 0=respect the periodic-data-only rule for final output derivatives (default), 1=post-blend FFT filtering
    int enable_initial_c3_effect_diagnostic = 1;              // 0=skip the first-step C3 effect CSV/print diagnostic, 1=preserve current behavior (default)
    int enable_nested_initialization_validation = 1;          // 0=skip the nested initialization validation bundle, 1=preserve current behavior (default)
};

inline constexpr int ActiveNestedL1Level = 1;

inline bool hasActiveNestedL1(
    const AlgorithmConfig& nested_config,
    int finest_level)
{
    return nested_config.use_crescent_algorithm == 1 &&
           finest_level >= ActiveNestedL1Level;
}

inline bool isActiveNestedL1Level(
    const AlgorithmConfig& nested_config,
    int lev)
{
    return nested_config.use_crescent_algorithm == 1 &&
           lev == ActiveNestedL1Level;
}

inline int resolveActiveBlendingMode(const AlgorithmConfig& nested_config)
{
    if (nested_config.blending_function_mode == 0 &&
        nested_config.use_sequential_blending == 1) {
        return 1;
    }

    return nested_config.blending_function_mode;
}

inline const char* blendingModeName(int blending_mode)
{
    switch (blending_mode) {
    case 0:
        return "Euclidean";
    case 1:
        return "sequential";
    case 2:
        return "harmonic";
    case 3:
        return "tensor-product";
    case 4:
        return "harmonic+septic";
    case 5:
        return "superellipse";
    default:
        return "unknown";
    }
}

} // namespace crescent::nested

#endif // CRESCENT_NESTED_CONFIG_H

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
    int state_interpolation_order = 2;              // L0->L1 eta/phi interpolation: 1=bilinear, 2=quadratic (default)
    int derivative_interpolation_order = 2;         // L0->L1 derivative interpolation: 1=bilinear, 2=quadratic (default)

    // C3 interface repair of the PTF periodic copies (CCF does not use it).
    amrex::Real c3_boundary_penalty = 0.01;           // endpoint penalty lambda_b
    amrex::Real c3_coefficient_regularization = 10.0; // coefficient regularization lambda_c

    int l1_coupling_mode = 1;                                 // 1=periodized total-field, PTF (default since 2026-08-20): tendency blending on C3-periodized fields. 2=compact correction-field, CCF: interpolated-parent background plus a tapered correction confined to the child modes above the parent's retained cutoff (docs/R2_CORRECTION_FIELD_SPEC.md). (the historical derivative blend, mode 0, and the unrestricted correction, ccf_band_restrict = 0, were retired in 2026-10)
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

} // namespace crescent::nested

#endif // CRESCENT_NESTED_CONFIG_H

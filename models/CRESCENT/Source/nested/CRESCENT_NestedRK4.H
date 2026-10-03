#ifndef CRESCENT_NESTED_RK4_H
#define CRESCENT_NESTED_RK4_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

/**
 * @file CRESCENT_NestedRK4.H
 * @brief Nested CRESCENT RK4 derivative blending and RHS utilities
 *
 * Implements the L1 derivative interpolation, core/L1 blending, and nonlinear
 * free-surface RHS pieces used by the nested CRESCENT RK4 stage driver.
 */

namespace crescent::nested::rk4 {

/**
 * @brief Interpolate core derivatives to the L1 grid at an RK4 substage time
 *
 * Performs both spatial interpolation (L0 grid → L1 grid) and temporal
 * interpolation (t^n and t^{n+1} → t^{n+c·Δt}) for all derivative fields.
 *
 * This function handles spatial + TEMPORAL interpolation for RK4 derivatives:
 *   D(t^n + c·Δt) = (1-c)·D^n + c·D^{n+1}
 *
 * Spatial interpolation order is selectable:
 *   1 = bilinear
 *   2 = quadratic
 *
 * Related function:
 *   - crescent::interpolation::interpolateParentToChild() handles spatial-only
 *     interpolation for state variables at initialization.
 *
 * Why separate?
 *   - This function: Spatial + temporal (needs two time levels for RK4)
 *   - Initialization helper: Spatial only (single time level)
 *
 * @param eta_x_L0_n    L0 ∂ₓη at time t^n
 * @param eta_x_L0_np1  L0 ∂ₓη at time t^{n+1}
 * @param c             RK4 substage coefficient (0, 0.5, 0.5, 1.0)
 * @param geom_L0       L0 geometry
 * @param geom_L1       L1 geometry
 * @param eta_x_interp  Output: interpolated ∂ₓη on L1 grid at substage time
 * @param deriv_stage_scratch Reused core-grid scratch for the temporal blend
 * @param coarse_patch_scratch Reused coarse-patch scratch for prolongation
 * @param interpolation_order Spatial interpolation order (1=bilinear, 2=quadratic)
 *
 * Note: This function handles ONE derivative component. Call for each of:
 *       (∂ₓη, ∂ᵧη, ∂ₓψ, ∂ᵧψ, w)
 */
void interpolateCoreDerivativeToL1Stage(
    const amrex::MultiFab& deriv_L0_n,
    const amrex::MultiFab& deriv_L0_np1,
    amrex::Real c,
    const amrex::Geometry& geom_L0,
    const amrex::Geometry& geom_L1,
    amrex::MultiFab& deriv_interp,
    amrex::MultiFab& deriv_stage_scratch,
    amrex::MultiFab& coarse_patch_scratch,
    int interpolation_order = 2
);

/**
 * @brief Blend a native L1 derivative with the core derivative interpolated to L1
 *
 * Blends the native L1 spectral derivative with the core derivative
 * interpolated to the L1 mesh:
 *   blended = (1-σ) · native_l1 + σ · core_on_l1
 *
 * Result:
 *   σ=0 in Ω_A → blended = native_l1  (pure L1 spectral)
 *   σ=1 in Ω_E → blended = core_on_l1 (pure L0 interpolation)
 *   σ∈(0,1) in Ω_M → smooth blend
 *
 * @param native_l1_derivative L1 spectral derivative on prepared periodic state
 * @param core_derivative_on_l1 L0 derivative interpolated to the L1 mesh
 * @param sigma     Blending function
 * @param blended_derivative Output: blended derivative used by the L1 RHS
 * @param correction_y_smoothing_passes Number of compact 1-2-1 smoothing passes
 *                                      applied to the correction field in y
 *                                      inside the mixing zone only
 */
void blendL1DerivativeWithCore(
    const amrex::MultiFab& native_l1_derivative,
    const amrex::MultiFab& core_derivative_on_l1,
    const amrex::MultiFab& sigma,
    amrex::MultiFab& blended_derivative,
    amrex::MultiFab& correction_delta_raw,
    amrex::MultiFab& correction_delta_smooth,
    int correction_y_smoothing_passes = 0
);

/**
 * @brief Step 5: Compute L1 RHS from the blended derivative bundle
 *
 * Computes RHS for L1 using BLENDED derivatives (from Step 4).
 * Full nonlinear equations consistent with L0 computeSurfaceRHS:
 *   ∂ₜη = -(∇ψ·∇η) + (1 + |∇η|²) · w
 *   ∂ₜψ = -g·η - ½|∇ψ|² + ½(1 + |∇η|²) · w²
 *
 * No intermediate FFT dealiasing is applied here. The blended L1 derivatives
 * are not periodic, so FFT-based dealiasing/filtering would violate the
 * periodic-only filtering rule used by CRESCENT.
 *
 * @param rhs_eta       Output: ∂ₜη (writes to component rk_stage)
 * @param rhs_phi       Output: ∂ₜψ (writes to component rk_stage)
 * @param eta           Current L1 surface elevation
 * @param eta_x         Blended ∂ₓη
 * @param eta_y         Blended ∂ᵧη
 * @param phi_x         Blended ∂ₓψ
 * @param phi_y         Blended ∂ᵧψ
 * @param phi_z         Blended vertical velocity
 * @param rk_stage      RK4 stage index (0-3)
 * @param gravity       Gravitational acceleration (default: 9.81)
 */
void computeSurfaceRHS_L1(
    amrex::MultiFab& rhs_eta,
    amrex::MultiFab& rhs_phi,
    const amrex::MultiFab& eta,
    const amrex::MultiFab& eta_x,
    const amrex::MultiFab& eta_y,
    const amrex::MultiFab& phi_x,
    const amrex::MultiFab& phi_y,
    const amrex::MultiFab& phi_z,
    int rk_stage,
    amrex::Real gravity = 9.81
);

} // namespace crescent::nested::rk4

#endif // CRESCENT_NESTED_RK4_H

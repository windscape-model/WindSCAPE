#ifndef CRESCENT_NESTED_C3PERIODICITY_H_
#define CRESCENT_NESTED_C3PERIODICITY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include "CRESCENT_NestedDomains.H"
#include "CRESCENT_NestedTypes.H"

/**
 * @file CRESCENT_NestedC3Periodicity.H
 * @brief C³-periodic seam repair for spectral methods on nested AMR
 * 
 * Instead of aggressively damping fields to zero at boundaries (which creates
 * spurious gradients), this approach makes fields C³-periodic by matching
 * value + first 3 derivatives at domain seams.
 * 
 * Key advantages:
 * - ROI (region of interest) completely untouched
 * - Minimal modification (only what's needed for periodicity)
 * - No spurious ∂chi_D/∂x or ∂chi_D/∂y terms in derivatives
 * - C∞ smooth at ROI/buffer interface
 */

namespace crescent::nested::c3 {

inline constexpr int InterfaceBasisOrder = 4;

/**
 * @brief Apply C³-periodic seam repair to make field spectral-compatible
 * 
 * Modifies field ONLY in buffer zones to enforce C³-periodicity:
 *   Δ₀ = u(L,y) - u(0,y) = 0          [value match]
 *   Δ₁ = u'(L,y) - u'(0,y) = 0        [1st derivative match]
 *   Δ₂ = u''(L,y) - u''(0,y) = 0      [2nd derivative match]
 *   Δ₃ = u'''(L,y) - u'''(0,y) = 0    [3rd derivative match]
 * 
 * ROI is preserved exactly (correction = 0 there).
 * 
 * @param field      MultiFab to make periodic (modified in place)
 * @param geom       Geometry for the level
 * @param ba         BoxArray for the level
 * @param buffer_width  Buffer zone width in grid points (recommend >= 16)
 * @param verbosity  Verbosity level for diagnostic output
 */
void applyC3SeamRepair(
    amrex::MultiFab& field,
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    int buffer_width = 16,
    int verbosity = 1
);

/**
 * @brief Apply C³-periodic seam repair using CRESCENT domain structure
 *
 * Uses the spectral domain boundaries from crescent_domains to determine
 * the C³-repair buffer zone. The repair is applied in the region between
 * the L1 box boundary and the spectral domain boundary (Ω_S).
 *
 * Ensures buffer width is at least min_buffer_width grid points.
 *
 * @param field      MultiFab to make periodic (modified in place)
 * @param geom       Geometry for the level
 * @param ba         BoxArray for the level
 * @param domains    CRESCENT domain structure with spectral boundaries
 * @param min_buffer_width  Minimum buffer width in grid points (default: 16)
 * @param verbosity  Verbosity level for diagnostic output
 */
void applyC3SeamRepair(
    amrex::MultiFab& field,
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    const crescent::nested::domains::DomainStructure& domains,
    int min_buffer_width = 16,
    int verbosity = 1
);

/**
 * @brief ACTIVE: Apply C³-periodic seam repair (interface-based, operates on extracted Ω_S)
 *
 * PRODUCTION VERSION - Used in CRESCENT nested AMR workflow.
 *
 * Algorithm:
 * - Operates on extracted Ω_S subdomain only
 * - Uses interface-based 5-point derivative reconstruction
 * - Polynomial basis Q_r(ξ) = (1-ξ)^4 · ξ^r embeds C³ continuity
 * - No separate rho mask (masking in polynomial structure)
 *
 * Workflow:
 *   1. Extract Ω_S from L1
 *   2. Apply C³-repair to Ω_S boundaries (makes Ω_S periodic)
 *   3. Insert C³-repaired Ω_S back into L1
 *
 * This is the validated production approach with excellent energy conservation.
 *
 * @param field      MultiFab on Ω_S subdomain (will be modified)
 * @param geom       Geometry for L1 (used for cell spacing)
 * @param ba         BoxArray for the field (Ω_S or L1)
 * @param buffer_width  Width of repair zones at Ω_S boundaries (in grid points)
 * @param verbosity  Verbosity level for diagnostic output
 * @param apply_x    Apply C³-repair in x-direction (default: true)
 * @param apply_y    Apply C³-repair in y-direction (default: true)
 */
struct C3InterfaceRepairOptions {
    int buffer_width = 12;
    int verbosity = 1;
    bool apply_x = true;
    bool apply_y = true;
    amrex::Real boundary_penalty_weight = 0.01;
    amrex::Real coefficient_regularization_weight = 1.0;
    int x_pass_y_smoothing_passes = 0;
    int x_pass_corner_taper = 0;
    int x_pass_alpha_y_smoothing_passes = 0;
    amrex::Real x_alpha_y_regularization_weight = 0.0;
    amrex::Real y_alpha_x_regularization_weight = 0.0;
    amrex::Real x_alpha_y_mean_regularization_weight = 0.0;
    amrex::Real y_alpha_x_mean_regularization_weight = 0.0;
    int x_transverse_harmonics = 0;
    int y_transverse_harmonics = 0;
    int x_alpha_y_projection_harmonics = 0;
    int y_alpha_x_projection_harmonics = 0;
};

/**
 * @brief Precompute the C3 interface basis (Qx, Qy, Bx, By) into a workspace.
 *
 * The basis depends only on the deterministic grid parameters
 * (InterfaceBasisOrder, nx, ny, buffer widths bx/by, and cell sizes dx/dy).
 * Call once at level setup; applyC3InterfaceSeamRepair then reads the cached
 * basis and does no per-call basis construction.
 */
void computeC3InterfaceBasis(
    crescent::nested::C3RepairWorkspace& workspace,
    int k,
    int nx,
    int ny,
    int bx,
    int by,
    amrex::Real dx,
    amrex::Real dy);

void applyC3InterfaceSeamRepair(
    amrex::MultiFab& field,
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    const C3InterfaceRepairOptions& options,
    crescent::nested::C3RepairWorkspace& workspace
);

void applyC3InterfaceSeamRepair(
    amrex::MultiFab& field,
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    int buffer_width,
    int verbosity,
    bool apply_x,
    bool apply_y,
    amrex::Real boundary_penalty_weight,
    amrex::Real coefficient_regularization_weight,
    int x_pass_y_smoothing_passes,
    int x_pass_corner_taper,
    int x_pass_alpha_y_smoothing_passes,
    amrex::Real x_alpha_y_regularization_weight,
    amrex::Real y_alpha_x_regularization_weight,
    amrex::Real x_alpha_y_mean_regularization_weight,
    amrex::Real y_alpha_x_mean_regularization_weight,
    int x_transverse_harmonics,
    int y_transverse_harmonics,
    int x_alpha_y_projection_harmonics,
    int y_alpha_x_projection_harmonics,
    crescent::nested::C3RepairWorkspace& workspace
);

/**
 * @brief Compute C∞ smooth window function (flat-top)
 *
 * Creates a buffer indicator: chi = 0 in ROI, chi = 1 in buffer zones,
 * with C∞ smooth transitions (all derivatives zero at transition points).
 *
 * @param chi_x      Output: 1D window in x-direction
 * @param chi_y      Output: 1D window in y-direction
 * @param nx, ny     Grid dimensions
 * @param b_x, b_y   Buffer width in grid points
 */
void computeSmoothWindow(
    amrex::Vector<amrex::Real>& chi_x,
    amrex::Vector<amrex::Real>& chi_y,
    int nx, int ny,
    int b_x, int b_y
);

} // namespace crescent::nested::c3

#endif // CRESCENT_NESTED_C3PERIODICITY_H_

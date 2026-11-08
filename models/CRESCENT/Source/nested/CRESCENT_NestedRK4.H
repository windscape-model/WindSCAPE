#ifndef CRESCENT_NESTED_RK4_H
#define CRESCENT_NESTED_RK4_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_LevelWaveOperatorView.H"

/**
 * @file CRESCENT_NestedRK4.H
 * @brief Tendency-blend (PTF) pieces of the child RK4 stage.
 */

namespace crescent::nested::rk4 {

/**
 * @brief blended = (1 - sigma) native_l1 + sigma parent_on_l1, pointwise.
 */
void blendWithParent(
    const amrex::MultiFab& native_l1_field,
    const amrex::MultiFab& parent_field_on_l1,
    const amrex::MultiFab& sigma,
    amrex::MultiFab& blended_field);

/**
 * @brief Native child RHS of the tendency blend
 *
 * Computes the full nonlinear free-surface RHS entirely from the C3-periodized
 * L1 fields (the repaired+filtered eta/phi copies, their spectral derivatives,
 * and the DNO vertical velocity), using the exact core recipe including 2/3
 * dealiasing of the nonlinear products. Every spectral operation acts on a
 * periodic field. Writes component 0 of the 1-component outputs.
 */
void computeNativeSurfaceRHS_L1(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const amrex::MultiFab& eta_periodic,
    const crescent::shared::DerivativeBundleView& native_bundle,
    const amrex::MultiFab& pressure,
    amrex::MultiFab& rhs_eta_native,
    amrex::MultiFab& rhs_phi_native
);

/**
 * @brief Tendency-blend coupling: interpolate L0's stage RHS to the L1 grid
 *
 * L0 advances before L1 within the same step, so component irk of the L0 RHS
 * MultiFab holds the exact L0 tendency at the RK4 substage time, so no
 * temporal interpolation is needed.
 */
void interpolateCoreStageRHSToL1(
    const amrex::MultiFab& rhs_core,
    int irk,
    const amrex::Geometry& geom_L0,
    const amrex::Geometry& geom_L1,
    const amrex::IntVect& ref_ratio,
    int interpolation_order,
    amrex::MultiFab& coarse_patch_scratch,
    amrex::MultiFab& rhs_core_on_l1
);

} // namespace crescent::nested::rk4

#endif // CRESCENT_NESTED_RK4_H

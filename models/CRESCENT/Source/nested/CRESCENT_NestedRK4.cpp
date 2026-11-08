#include "CRESCENT_NestedRK4.H"

// Tendency-blend (PTF) helpers of the child RK4 stage.
#include "CRESCENT_Constants.H"
#include "CRESCENT_Dynamics.H"
#include "CRESCENT_Interpolation.H"

namespace crescent::nested::rk4 {

using namespace amrex;

void blendWithParent(
    const MultiFab& native_l1_field,
    const MultiFab& parent_field_on_l1,
    const MultiFab& sigma,
    MultiFab& blended_field)
{
    // B[native_l1, parent_on_l1] = (1-sigma) native_l1 + sigma parent_on_l1
    for (MFIter mfi(blended_field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto native_arr = native_l1_field.const_array(mfi);
        auto core_arr = parent_field_on_l1.const_array(mfi);
        auto sigma_arr = sigma.const_array(mfi);
        auto blended_arr = blended_field.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real s = sigma_arr(i,j,k);
            blended_arr(i,j,k) = (1.0 - s) * native_arr(i,j,k) + s * core_arr(i,j,k);
        });
    }
}

void computeNativeSurfaceRHS_L1(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const MultiFab& eta_periodic,
    const crescent::shared::DerivativeBundleView& native_bundle,
    const MultiFab& pressure,
    MultiFab& rhs_eta_native,
    MultiFab& rhs_phi_native)
{
    // The core RHS kernel is reusable on L1 because its spectral pieces
    // (dealiasing of metric, rhs_eta, phi_z^2, rhs_phi) derive their grid
    // metrics from the MultiFab minimalBox and the level cell size, both of
    // which are the L1 patch values here; only the FFT object must be the
    // L1 one. All inputs are the C3-periodized fields, so the 2/3 truncation
    // is exact for the quadratic products, as on L0.
    const int lev = level_wave_operator.lev;
    crescent::dynamics::computeSurfaceRHS(
        rhs_eta_native,
        rhs_phi_native,
        eta_periodic,
        native_bundle.eta_x,
        native_bundle.eta_y,
        native_bundle.phi_x,
        native_bundle.phi_y,
        native_bundle.phi_z,
        pressure,
        level_wave_operator.spectral_workspace.compute_rhs_metric[lev],
        level_wave_operator.spectral_workspace.compute_rhs_phi_z_sq[lev],
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        crescent::constants::INV_GRAVITY,
        level_wave_operator.perturbation_order,
        0,
        level_wave_operator.spectral_workspace.dealiasing_scratch[lev]);
}

void interpolateCoreStageRHSToL1(
    const MultiFab& rhs_core,
    int irk,
    const Geometry& geom_L0,
    const Geometry& geom_L1,
    const IntVect& ref_ratio,
    int interpolation_order,
    MultiFab& coarse_patch_scratch,
    MultiFab& rhs_core_on_l1)
{
    MultiFab rhs_core_stage(
        const_cast<MultiFab&>(rhs_core), amrex::make_alias, irk, 1);

    crescent::interpolation::interpolateScalarParentToChild(
        rhs_core_stage,
        rhs_core_on_l1,
        geom_L0,
        geom_L1,
        ref_ratio,
        coarse_patch_scratch,
        interpolation_order);
}

} // namespace crescent::nested::rk4

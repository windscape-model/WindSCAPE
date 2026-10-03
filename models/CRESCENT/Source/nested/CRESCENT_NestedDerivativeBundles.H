#ifndef CRESCENT_NESTED_DERIVATIVE_BUNDLES_H
#define CRESCENT_NESTED_DERIVATIVE_BUNDLES_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>

#include "CRESCENT_NestedRK4.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_NestedTypes.H"

namespace crescent::nested::bundles {

inline StageDerivativeViews makeL1StageDerivativeViews(
    CouplingDerivativeState& derivative_coupling,
    int lev)
{
    return {
        derivative_coupling.core_snapshots.n.bundleAt(0),
        derivative_coupling.core_snapshots.np1.bundleAt(0),
        derivative_coupling.core_derivatives_on_l1.bundleAt(lev),
        derivative_coupling.l1_blended_derivatives.bundleAt(lev)
    };
}

inline OutputDerivativeViews makeL1OutputDerivativeViews(
    CouplingDerivativeState& derivative_coupling,
    int lev)
{
    return {
        derivative_coupling.core_snapshots.np1.bundleAt(0),
        derivative_coupling.core_derivatives_on_l1.bundleAt(lev),
        derivative_coupling.l1_blended_derivatives.bundleAt(lev)
    };
}

inline LinearDerivativeSources selectLinearDerivativeSources(
    const PeriodicDerivativeView& periodic_fields,
    bool use_filtered_periodic_copies)
{
    return use_filtered_periodic_copies ?
        periodic_fields.periodicState() :
        periodic_fields.linearPrefilterState();
}

inline void copyDerivativeBundle(shared::DerivativeBundleView dst_bundle,
                                 shared::DerivativeBundleView src_bundle)
{
    amrex::MultiFab::Copy(dst_bundle.eta_x, src_bundle.eta_x, 0, 0, 1, 0);
    amrex::MultiFab::Copy(dst_bundle.eta_y, src_bundle.eta_y, 0, 0, 1, 0);
    amrex::MultiFab::Copy(dst_bundle.phi_x, src_bundle.phi_x, 0, 0, 1, 0);
    amrex::MultiFab::Copy(dst_bundle.phi_y, src_bundle.phi_y, 0, 0, 1, 0);
    amrex::MultiFab::Copy(dst_bundle.phi_z, src_bundle.phi_z, 0, 0, 1, 0);
}

inline void fillDerivativeBundleBoundaries(shared::DerivativeBundleView bundle,
                                           const amrex::Periodicity& periodicity)
{
    bundle.eta_x.FillBoundary(periodicity);
    bundle.eta_y.FillBoundary(periodicity);
    bundle.phi_x.FillBoundary(periodicity);
    bundle.phi_y.FillBoundary(periodicity);
    bundle.phi_z.FillBoundary(periodicity);
}

inline void interpolateCoreDerivativesToL1(
    shared::DerivativeBundleView core_bundle_n,
    shared::DerivativeBundleView core_bundle_np1,
    amrex::Real temporal_weight,
    const amrex::Geometry& geom_core,
    const amrex::Geometry& geom_nested,
    int interpolation_order,
    amrex::MultiFab& core_stage_derivative_scratch,
    amrex::MultiFab& core_interp_patch_scratch,
    shared::DerivativeBundleView core_derivatives_on_l1_bundle)
{
    rk4::interpolateCoreDerivativeToL1Stage(
        core_bundle_n.eta_x, core_bundle_np1.eta_x, temporal_weight,
        geom_core, geom_nested, core_derivatives_on_l1_bundle.eta_x,
        core_stage_derivative_scratch, core_interp_patch_scratch,
        interpolation_order);
    rk4::interpolateCoreDerivativeToL1Stage(
        core_bundle_n.eta_y, core_bundle_np1.eta_y, temporal_weight,
        geom_core, geom_nested, core_derivatives_on_l1_bundle.eta_y,
        core_stage_derivative_scratch, core_interp_patch_scratch,
        interpolation_order);
    rk4::interpolateCoreDerivativeToL1Stage(
        core_bundle_n.phi_x, core_bundle_np1.phi_x, temporal_weight,
        geom_core, geom_nested, core_derivatives_on_l1_bundle.phi_x,
        core_stage_derivative_scratch, core_interp_patch_scratch,
        interpolation_order);
    rk4::interpolateCoreDerivativeToL1Stage(
        core_bundle_n.phi_y, core_bundle_np1.phi_y, temporal_weight,
        geom_core, geom_nested, core_derivatives_on_l1_bundle.phi_y,
        core_stage_derivative_scratch, core_interp_patch_scratch,
        interpolation_order);
    rk4::interpolateCoreDerivativeToL1Stage(
        core_bundle_n.phi_z, core_bundle_np1.phi_z, temporal_weight,
        geom_core, geom_nested, core_derivatives_on_l1_bundle.phi_z,
        core_stage_derivative_scratch, core_interp_patch_scratch,
        interpolation_order);
}

inline void blendL1DerivativesWithCore(shared::DerivativeBundleView native_l1_bundle,
                                       shared::DerivativeBundleView core_derivatives_on_l1_bundle,
                                       const amrex::MultiFab& sigma_blend,
                                       int correction_y_smoothing_passes,
                                       amrex::MultiFab& correction_delta_raw,
                                       amrex::MultiFab& correction_delta_smooth,
                                       shared::DerivativeBundleView blended_derivative_bundle)
{
    rk4::blendL1DerivativeWithCore(
        native_l1_bundle.eta_x, core_derivatives_on_l1_bundle.eta_x, sigma_blend,
        blended_derivative_bundle.eta_x,
        correction_delta_raw,
        correction_delta_smooth,
        correction_y_smoothing_passes);
    rk4::blendL1DerivativeWithCore(
        native_l1_bundle.eta_y, core_derivatives_on_l1_bundle.eta_y, sigma_blend,
        blended_derivative_bundle.eta_y,
        correction_delta_raw,
        correction_delta_smooth,
        correction_y_smoothing_passes);
    rk4::blendL1DerivativeWithCore(
        native_l1_bundle.phi_x, core_derivatives_on_l1_bundle.phi_x, sigma_blend,
        blended_derivative_bundle.phi_x,
        correction_delta_raw,
        correction_delta_smooth,
        correction_y_smoothing_passes);
    rk4::blendL1DerivativeWithCore(
        native_l1_bundle.phi_y, core_derivatives_on_l1_bundle.phi_y, sigma_blend,
        blended_derivative_bundle.phi_y,
        correction_delta_raw,
        correction_delta_smooth,
        correction_y_smoothing_passes);
    rk4::blendL1DerivativeWithCore(
        native_l1_bundle.phi_z, core_derivatives_on_l1_bundle.phi_z, sigma_blend,
        blended_derivative_bundle.phi_z,
        correction_delta_raw,
        correction_delta_smooth,
        correction_y_smoothing_passes);
}

} // namespace crescent::nested::bundles

#endif // CRESCENT_NESTED_DERIVATIVE_BUNDLES_H

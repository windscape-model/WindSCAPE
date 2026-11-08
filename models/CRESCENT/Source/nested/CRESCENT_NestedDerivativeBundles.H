#ifndef CRESCENT_NESTED_DERIVATIVE_BUNDLES_H
#define CRESCENT_NESTED_DERIVATIVE_BUNDLES_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Periodicity.H>

#include "CRESCENT_Interpolation.H"
#include "CRESCENT_NestedRK4.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_NestedTypes.H"

namespace crescent::nested::bundles {

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

// Prolong each field of a parent bundle to the child grid and fill the
// child ghost cells periodically.
inline void interpolateParentBundleToChild(
    shared::DerivativeBundleView parent_bundle,
    const amrex::Geometry& parent_geom,
    const amrex::Geometry& child_geom,
    const amrex::IntVect& ref_ratio,
    int interpolation_order,
    amrex::MultiFab& coarse_patch_scratch,
    shared::DerivativeBundleView child_bundle)
{
    amrex::MultiFab* const parent_fields[] = {
        &parent_bundle.eta_x, &parent_bundle.eta_y, &parent_bundle.phi_x,
        &parent_bundle.phi_y, &parent_bundle.phi_z};
    amrex::MultiFab* const child_fields[] = {
        &child_bundle.eta_x, &child_bundle.eta_y, &child_bundle.phi_x,
        &child_bundle.phi_y, &child_bundle.phi_z};
    for (int f = 0; f < 5; ++f) {
        crescent::interpolation::interpolateScalarParentToChild(
            *parent_fields[f], *child_fields[f], parent_geom, child_geom,
            ref_ratio, coarse_patch_scratch, interpolation_order);
        child_fields[f]->FillBoundary(child_geom.periodicity());
    }
}

// blended = (1 - sigma) native + sigma parent_on_child, field by field.
inline void blendBundleWithParent(
    shared::DerivativeBundleView native_l1_bundle,
    shared::DerivativeBundleView parent_on_l1_bundle,
    const amrex::MultiFab& sigma_blend,
    shared::DerivativeBundleView blended_bundle)
{
    rk4::blendWithParent(native_l1_bundle.eta_x, parent_on_l1_bundle.eta_x,
                         sigma_blend, blended_bundle.eta_x);
    rk4::blendWithParent(native_l1_bundle.eta_y, parent_on_l1_bundle.eta_y,
                         sigma_blend, blended_bundle.eta_y);
    rk4::blendWithParent(native_l1_bundle.phi_x, parent_on_l1_bundle.phi_x,
                         sigma_blend, blended_bundle.phi_x);
    rk4::blendWithParent(native_l1_bundle.phi_y, parent_on_l1_bundle.phi_y,
                         sigma_blend, blended_bundle.phi_y);
    rk4::blendWithParent(native_l1_bundle.phi_z, parent_on_l1_bundle.phi_z,
                         sigma_blend, blended_bundle.phi_z);
}

} // namespace crescent::nested::bundles

#endif // CRESCENT_NESTED_DERIVATIVE_BUNDLES_H

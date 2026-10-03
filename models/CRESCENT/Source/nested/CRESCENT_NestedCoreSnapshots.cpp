#include "CRESCENT_NestedCoreSnapshots.H"

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedDerivativeBundles.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"

namespace {

void computeAndFillCoreDerivativeBundle(
    int lev,
    const crescent::ownership::CoreSolverAccess& core,
    const amrex::Vector<amrex::Geometry>& geom,
    crescent::shared::DerivativeBundleView core_bundle_out)
{
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            geom[lev].data());
    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        core,
        core_bundle_out);

    crescent::nested::bundles::fillDerivativeBundleBoundaries(
        core_bundle_out,
        geom[lev].periodicity());
}

} // namespace

namespace crescent::nested::snapshots {

void updateCoreInterpolationSourceAtStepStart(
    int finest_level,
    int lev,
    const crescent::ownership::CoreSolverAccess& core,
    const crescent::ownership::NestedSolverAccess& nested,
    const amrex::Vector<amrex::Geometry>& geom)
{
    if (!crescent::nested::hasActiveNestedL1(nested.config, finest_level) ||
        lev != 0) {
        return;
    }

    auto& coupling_derivatives = nested.runtime.coupling_derivatives;
    auto core_interp_source_n_bundle = coupling_derivatives.core_snapshots.n.bundleAt(0);
    auto core_interp_source_np1_bundle = coupling_derivatives.core_snapshots.np1.bundleAt(0);

    if (core.runtime.istep[0] > 0) {
        crescent::nested::bundles::copyDerivativeBundle(
            core_interp_source_n_bundle,
            core_interp_source_np1_bundle);
        return;
    }

    computeAndFillCoreDerivativeBundle(
        0,
        core,
        geom,
        core_interp_source_n_bundle);
    crescent::nested::bundles::copyDerivativeBundle(
        core_interp_source_np1_bundle,
        core_interp_source_n_bundle);
    crescent::nested::bundles::fillDerivativeBundleBoundaries(
        core_interp_source_np1_bundle,
        geom[0].periodicity());
}

void storeCoreInterpolationSourceAtStepEnd(
    int finest_level,
    int lev,
    const crescent::ownership::CoreSolverAccess& core,
    const crescent::ownership::NestedSolverAccess& nested,
    const amrex::Vector<amrex::Geometry>& geom)
{
    if (lev != 0 ||
        !crescent::nested::hasActiveNestedL1(nested.config, finest_level)) {
        return;
    }

    auto& coupling_derivatives = nested.runtime.coupling_derivatives;
    auto core_final_bundle = coupling_derivatives.core_snapshots.np1.bundleAt(0);
    auto core_final_output_bundle = core.solution.derivativeBundleAt(0);

    computeAndFillCoreDerivativeBundle(
        lev,
        core,
        geom,
        core_final_bundle);
    crescent::nested::bundles::copyDerivativeBundle(
        core_final_output_bundle,
        core_final_bundle);
}

} // namespace crescent::nested::snapshots

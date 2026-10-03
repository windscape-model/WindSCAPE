#include "CRESCENT_RestartRebuild.H"

#include <AMReX_IntVect.H>

#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedCoreSnapshots.H"
#include "CRESCENT_NestedDerivativeBundles.H"
#include "CRESCENT_NestedFreeSurfaceEvolution.H"
#include "CRESCENT_NestedSetupDerivedFields.H"
#include "CRESCENT_NestedSetupStorage.H"
#include "CRESCENT_PostStep.H"
#include "CRESCENT_SetupDerivedFields.H"
#include "CRESCENT_SetupLifecycle.H"
#include "CRESCENT_SetupSpectralWorkspace.H"
#include "CRESCENT_SetupWaveState.H"
#include "CRESCENT_SolverInitializationContext.H"

namespace crescent::solver::checkpoint {

void rebuildRestartDerivedState(
    const crescent::solver::setup::SolverInitializationContext& initialization_context)
{
    auto core = initialization_context.coreAccess();
    auto nested = initialization_context.nestedAccess();
    const auto mesh = initialization_context.meshHierarchy();
    const int finest_level = initialization_context.finestLevel();

    crescent::solver::setup::printRefinementSummary(
        core.runtime,
        mesh,
        finest_level);
    crescent::solver::setup::initializeSpectralWorkspaces(
        initialization_context);

    for (int lev = 0; lev <= finest_level; ++lev) {
        crescent::solver::setup::wave_state::initializeLevelSpectralMetrics(
            initialization_context.levelContext(lev));
        core.solution.eta[lev].FillBoundary(mesh.geom[lev].periodicity());
        core.solution.phi[lev].FillBoundary(mesh.geom[lev].periodicity());
        core.solution.pressure[lev].FillBoundary(mesh.geom[lev].periodicity());
    }

    crescent::solver::setup::derived_fields::initializeCoreDerivedFields(
        initialization_context.levelContext(0));

    if (!crescent::nested::hasActiveNestedL1(nested.config, finest_level)) {
        return;
    }

    const int lev = crescent::nested::ActiveNestedL1Level;
    const auto l1_initialization_context = initialization_context.levelContext(lev);
    crescent::nested::setup::storage::initializeDomainsAndStorage(
        l1_initialization_context);
    crescent::nested::setup::derived_fields::initializeFromState(
        l1_initialization_context);

    crescent::nested::snapshots::storeCoreInterpolationSourceAtStepEnd(
        finest_level,
        0,
        core,
        nested,
        mesh.geom);
    crescent::nested::bundles::copyDerivativeBundle(
        nested.runtime.coupling_derivatives.core_snapshots.n.bundleAt(0),
        nested.runtime.coupling_derivatives.core_snapshots.np1.bundleAt(0));

    const amrex::IntVect core_to_l1_ref_ratio =
        mesh.refinementRatio(lev);
    auto l1_evolution_context = crescent::shared::makeLevelEvolutionContext(
        lev,
        finest_level,
        core,
        nested,
        mesh,
        core_to_l1_ref_ratio);
    crescent::nested::free_surface_evolution::recomputeL1OutputDerivatives(
        l1_evolution_context);
    crescent::solver::poststep::computeVelocityOutputs(
        l1_evolution_context);
}

} // namespace crescent::solver::checkpoint

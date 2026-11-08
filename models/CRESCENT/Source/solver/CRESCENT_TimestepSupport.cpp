#include "CRESCENT_FreeSurfaceEvolution.H"

#include <functional>

#include "CRESCENT_CoreFreeSurfaceEvolution.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedCoreSnapshots.H"
#include "CRESCENT_NestedCorrectionField.H"
#include "CRESCENT_NestedTimestep.H"
#include "CRESCENT_PostStep.H"
#include "CRESCENT_Spectral.H"

namespace crescent::solver::free_surface_evolution {

namespace {


void advanceLevelTimestep(
    crescent::shared::LevelContext& level_context,
    amrex::Real time);

bool usesCorrectionFieldCoupling(const crescent::shared::LevelContext& level_context)
{
    return level_context.nestedConfig().l1_coupling_mode == 2 &&
           crescent::nested::hasActiveNestedL1(
               level_context.nestedConfig(),
               level_context.finestLevel());
}

void advanceActiveLevelTimestepPath(
    crescent::shared::LevelContext& level_context,
    amrex::Real time)
{
    const int lev = level_context.level();
    if (lev == 0) {
        if (usesCorrectionFieldCoupling(level_context)) {
            const std::function<void(int)> capture_hook =
                [&level_context](int rk_stage) {
                    crescent::nested::correction_field::captureCoreStage(
                        level_context,
                        rk_stage);
                };
            crescent::core::free_surface_evolution::advanceLevelTimestep(
                level_context,
                &capture_hook);
        } else {
            crescent::core::free_surface_evolution::advanceLevelTimestep(
                level_context);
        }
        return;
    }

    if (crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev)) {
        crescent::nested::free_surface_evolution::advanceL1Timestep(
            level_context,
            time);
        return;
    }

    amrex::Abort("Nested L1 runs require crescent.use_crescent_algorithm = 1");
}

void refreshCoreInterpolationSourcesAtStepEnd(
    crescent::shared::LevelContext& level_context)
{
    crescent::nested::snapshots::storeCoreInterpolationSourceAtStepEnd(
        level_context.finestLevel(),
        level_context.level(),
        level_context.core(),
        level_context.nested(),
        level_context.meshHierarchy().geom);
}

void finalizeSurfaceLevelOutputs(
    crescent::shared::LevelContext& level_context)
{
    // The standard derivative bundle is needed only on the parent of a
    // single-level run: with a child, the step-end parent snapshot has
    // already computed it from the same state and copied it into the
    // solution bundle (storeCoreInterpolationSourceAtStepEnd), and the child
    // output bundles come from the coupling (the live child state is
    // non-periodic and is never transformed).
    if (!crescent::nested::hasActiveNestedL1(
            level_context.nestedConfig(),
            level_context.finestLevel())) {
        crescent::solver::poststep::recomputeStandardDerivativeBundle(
            level_context);
    }

    crescent::solver::poststep::computeVelocityOutputs(
        level_context);
}

void advanceLevelTimestep(
    crescent::shared::LevelContext& level_context,
    amrex::Real time)
{
    advanceActiveLevelTimestepPath(
        level_context,
        time);

    refreshCoreInterpolationSourcesAtStepEnd(level_context);

    crescent::nested::free_surface_evolution::recomputeL1OutputDerivatives(
        level_context);

    finalizeSurfaceLevelOutputs(level_context);
}

} // namespace

void advanceAllLevelsForStep(
    crescent::shared::SolverContext& evolution_context,
    amrex::Real& time)
{
    for (int lev = 0; lev <= evolution_context.finestLevel(); ++lev) {
        auto level_context = evolution_context.levelContext(lev);
        advanceLevelTimestep(level_context, time);
    }
}

void updateTimeStateAfterStep(
    crescent::shared::SolverContext& evolution_context,
    amrex::Real& time)
{
    auto& runtime = evolution_context.runtime();
    for (int lev = 0; lev <= evolution_context.finestLevel(); ++lev) {
        runtime.t_new[lev] = time + runtime.dt[lev];
        runtime.istep[lev] += 1;
    }
    time += runtime.dt[0];
}

} // namespace crescent::solver::free_surface_evolution

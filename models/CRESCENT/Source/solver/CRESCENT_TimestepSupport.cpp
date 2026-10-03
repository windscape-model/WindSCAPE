#include "CRESCENT_FreeSurfaceEvolution.H"

#include "CRESCENT_CoreFreeSurfaceEvolution.H"
#include "CRESCENT_FreeSurfaceEvolutionContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedCoreSnapshots.H"
#include "CRESCENT_NestedFreeSurfaceEvolution.H"
#include "CRESCENT_PostStep.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_TimestepLogging.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"
#include "shared/CRESCENT_LevelWaveOperatorView.H"

namespace crescent::solver::free_surface_evolution {

namespace {

using FreeSurfaceEvolutionContext = crescent::solver::FreeSurfaceEvolutionContext;
using LevelEvolutionContext = crescent::shared::LevelEvolutionContext;

void advanceLevelTimestep(
    crescent::shared::LevelEvolutionContext& level_context,
    amrex::Real time);

void advanceActiveLevelTimestepPath(
    LevelEvolutionContext& level_context,
    amrex::Real time)
{
    const int lev = level_context.level();
    if (lev == 0) {
        crescent::core::free_surface_evolution::advanceLevelTimestep(
            level_context);
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

void logAdvanceStart(
    const LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    const auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();

    if (runtime.verbosity >= 2) {
        amrex::Print() << "\n=== ADVANCE_WAVE START (timestep "
                       << runtime.istep[lev] + 1
                       << ") ===" << std::endl;
    }
    crescent::solver::logging::logAdvanceWaveStateSnapshot(
        solution,
        runtime,
        lev,
        "Before RK4");
}

void refreshCoreInterpolationSourcesAtStepStart(
    LevelEvolutionContext& level_context)
{
    crescent::nested::snapshots::updateCoreInterpolationSourceAtStepStart(
        level_context.finestLevel(),
        level_context.level(),
        level_context.coreAccess(),
        level_context.nestedAccess(),
        level_context.meshHierarchy().geom);
}

void refreshCoreInterpolationSourcesAtStepEnd(
    LevelEvolutionContext& level_context)
{
    crescent::nested::snapshots::storeCoreInterpolationSourceAtStepEnd(
        level_context.finestLevel(),
        level_context.level(),
        level_context.coreAccess(),
        level_context.nestedAccess(),
        level_context.meshHierarchy().geom);
}

void finalizeSurfaceLevelOutputs(
    LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    const auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());
    const auto& spectral_settings = level_wave_operator.spectral_settings;

    crescent::solver::logging::logAdvanceWaveStateSnapshot(
        solution,
        runtime,
        lev,
        "\nAfter surf_update");
    crescent::solver::logging::logAdvanceWavePreFilterDiagnostics(
        solution,
        runtime,
        lev,
        spectral_settings.filter.beta1_filter,
        spectral_settings.filter.beta2_filter);

    crescent::solver::poststep::applyOptionalPostTimestepLowPassFilter(
        level_context);

    crescent::solver::logging::logAdvanceWavePostFilterDiagnostics(
        solution,
        runtime,
        lev);

    crescent::solver::poststep::recomputeStandardDerivativeBundle(
        level_context);

    crescent::solver::logging::logAdvanceWaveDerivativeRecomputeDiagnostics(
        solution,
        runtime,
        lev);

    crescent::solver::logging::logAdvanceWaveEndBanner(runtime);

    crescent::solver::poststep::computeVelocityOutputs(
        level_context);
    crescent::solver::logging::logFinalWaveTimestepState(
        solution,
        runtime,
        lev);
}

void advanceLevelTimestep(
    crescent::shared::LevelEvolutionContext& level_context,
    amrex::Real time)
{
    logAdvanceStart(level_context);
    refreshCoreInterpolationSourcesAtStepStart(level_context);

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
    crescent::solver::FreeSurfaceEvolutionContext& evolution_context,
    amrex::Real& time)
{
    for (int lev = 0; lev <= evolution_context.finestLevel(); ++lev) {
        auto level_context = evolution_context.levelContext(lev);
        advanceLevelTimestep(level_context, time);
    }
}

void updateTimeStateAfterStep(
    crescent::solver::FreeSurfaceEvolutionContext& evolution_context,
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

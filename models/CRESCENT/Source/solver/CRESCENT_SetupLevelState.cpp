#include "CRESCENT_SetupLevelState.H"

#include "CRESCENT_LevelInitializationContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedSetup.H"
#include "CRESCENT_LevelInitLogging.H"
#include "CRESCENT_SetupDerivedFields.H"
#include "CRESCENT_SetupWaveState.H"

namespace {

void finalizeCoreLevelInitialization(
    const crescent::shared::LevelInitializationContext& level_context)
{
    crescent::solver::setup::derived_fields::initializeCoreDerivedFields(
        level_context);
    level_context.solutionState().pressure[0].setVal(0.0);
}

void finalizeNestedLevelInitialization(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev > 0);
    auto& solution = level_context.solutionState();
    solution.pressure[lev].setVal(0.0);

    if (crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev)) {
        crescent::nested::setup::initializeNestedL1Level(
            level_context);
    }
}

void initializeLevelState(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    auto& runtime = level_context.runtime();
    auto& wave = level_context.wave();
    crescent::solver::setup::wave_state::prepareLevelWaveState(
        level_context);
    crescent::solver::logging::logLevelInitializationState(
        solution,
        runtime,
        wave,
        lev,
        "Before");
    crescent::solver::setup::wave_state::applyInitializationWaveFilter(
        level_context);
    crescent::solver::logging::logLevelInitializationState(
        solution,
        runtime,
        wave,
        lev,
        "After");
    if (lev == 0) {
        finalizeCoreLevelInitialization(level_context);
    } else {
        finalizeNestedLevelInitialization(level_context);
    }
}

} // namespace

namespace crescent::solver::setup {

void initializeLevelStates(
    const SolverInitializationContext& initialization_context)
{
    for (int lev = 0; lev <= initialization_context.finestLevel(); ++lev) {
        initializeLevelState(
            initialization_context.levelContext(lev));
    }
}

} // namespace crescent::solver::setup

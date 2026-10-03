#pragma once

#include "CRESCENT_LevelInitializationContext.H"

namespace crescent::solver::setup::wave_state {

void prepareLevelWaveState(
    const crescent::shared::LevelInitializationContext& level_context);

void initializeLevelSpectralMetrics(
    const crescent::shared::LevelInitializationContext& level_context);

void applyInitializationWaveFilter(
    const crescent::shared::LevelInitializationContext& level_context);

} // namespace crescent::solver::setup::wave_state

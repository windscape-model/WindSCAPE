#pragma once

#include "CRESCENT_LevelMeshView.H"
#include "CRESCENT_SolverState.H"

namespace crescent::init {

// Initialize the level-0 wave field and its wave-dependent runtime metadata.
void initializeCoreWaveField(
    CoreSolverState& core,
    shared::LevelMeshView level_mesh,
    int finest_level);

} // namespace crescent::init

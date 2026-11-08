#pragma once

#include <string>

#include "CRESCENT_LevelMeshView.H"

namespace crescent {
struct CoreSolverState;
struct SolverState;
}

namespace crescent::shared {
class SolverContext;
}

namespace crescent::solver::setup {

// Reads the inputs and sizes the per-level containers; called by the
// CRESCENT constructor before AMReX builds the level hierarchy. The solver
// inputs are read under solver_pp_prefix ("crescent" unless an in-process
// driver gives the instance its own prefix).
void constructSolverState(
    crescent::SolverState& state,
    int max_level,
    const std::string& solver_pp_prefix = "crescent");

// Defines the field storage of one level (zero-initialized); the wave
// initializers or the checkpoint reader fill eta and phi.
void initializeScratchLevelStorage(
    int lev,
    crescent::CoreSolverState& core,
    shared::LevelMeshView level_mesh);

// Fresh run, after AMReX built the hierarchy: spectral workspaces, the wave
// field on the parent, the child interpolated from it, derived fields.
void initializeSolverState(
    const crescent::shared::SolverContext& context);

// Restart, after the checkpoint reader restored eta, phi and pressure: the
// same workspaces and derived fields, rebuilt from the restored state.
void rebuildRestartDerivedState(
    const crescent::shared::SolverContext& context);

} // namespace crescent::solver::setup

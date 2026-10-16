#ifndef CRESCENT_SOLVER_STATE_H
#define CRESCENT_SOLVER_STATE_H

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_WaveState.H"

namespace crescent {

// Run control, wave parameters, per-level fields and spectral workspaces.
struct CoreSolverState
{
    shared::RuntimeControlState runtime;
    shared::WaveRuntimeState wave;
    shared::SolutionState solution;
    shared::RootLayouts root_layouts;
    shared::SpectralWorkspace spectral_workspace;
};

// Nested-coupling configuration and the child-level coupling state.
struct NestedSolverState
{
    nested::AlgorithmConfig config;
    nested::RuntimeState runtime;
};

// All solver state owned by the CRESCENT AmrCore object.
struct SolverState
{
    CoreSolverState core;
    NestedSolverState nested;
};

} // namespace crescent

#endif // CRESCENT_SOLVER_STATE_H

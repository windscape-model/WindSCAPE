#ifndef CRESCENT_SETUP_INITIALIZATION_H
#define CRESCENT_SETUP_INITIALIZATION_H

#include "CRESCENT_SolverInitializationContext.H"
#include "CRESCENT_SetupLifecycle.H"
#include "CRESCENT_SetupLevelState.H"
#include "CRESCENT_SetupSpectralWorkspace.H"

namespace crescent::solver::setup {

inline void initializeSolverState(
    const SolverInitializationContext& initialization_context)
{
    printRefinementSummary(
        initialization_context.runtime(),
        initialization_context.meshHierarchy(),
        initialization_context.finestLevel());
    initializeSpectralWorkspaces(initialization_context);
    initializeLevelStates(initialization_context);
}

} // namespace crescent::solver::setup

#endif // CRESCENT_SETUP_INITIALIZATION_H

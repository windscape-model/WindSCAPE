#ifndef CRESCENT_FREE_SURFACE_EVOLUTION_H
#define CRESCENT_FREE_SURFACE_EVOLUTION_H

#include <AMReX_REAL.H>

namespace crescent::shared {
class SolverContext;
}

namespace crescent::solver::free_surface_evolution {

void advanceAllLevelsForStep(
    crescent::shared::SolverContext& evolution_context,
    amrex::Real& time);

void updateTimeStateAfterStep(
    crescent::shared::SolverContext& evolution_context,
    amrex::Real& time);

} // namespace crescent::solver::free_surface_evolution

#endif // CRESCENT_FREE_SURFACE_EVOLUTION_H

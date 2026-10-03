#ifndef CRESCENT_NESTED_FREE_SURFACE_EVOLUTION_H
#define CRESCENT_NESTED_FREE_SURFACE_EVOLUTION_H

#include <AMReX_Real.H>

namespace crescent::shared {
class LevelEvolutionContext;
}

namespace crescent::nested::free_surface_evolution {

void advanceL1Timestep(
    crescent::shared::LevelEvolutionContext& level_context,
    amrex::Real time);

void recomputeL1OutputDerivatives(
    crescent::shared::LevelEvolutionContext& level_context);

} // namespace crescent::nested::free_surface_evolution

#endif // CRESCENT_NESTED_FREE_SURFACE_EVOLUTION_H

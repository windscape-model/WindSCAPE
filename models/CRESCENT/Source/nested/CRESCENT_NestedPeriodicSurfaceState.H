#ifndef CRESCENT_NESTED_PERIODIC_SURFACE_STATE_H
#define CRESCENT_NESTED_PERIODIC_SURFACE_STATE_H

#include "CRESCENT_NestedFreeSurfaceEvolution.H"

namespace crescent::nested::periodic_surface_state {

void applyPostRK4CoreStateBlend(
    crescent::shared::LevelEvolutionContext& level_context);

} // namespace crescent::nested::periodic_surface_state

#endif // CRESCENT_NESTED_PERIODIC_SURFACE_STATE_H

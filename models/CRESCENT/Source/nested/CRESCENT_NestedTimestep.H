#ifndef CRESCENT_NESTED_TIMESTEP_H
#define CRESCENT_NESTED_TIMESTEP_H

#include <AMReX_REAL.H>

namespace crescent::shared {
class LevelContext;
}

namespace crescent::nested::free_surface_evolution {

// One child step after the parent step: four RK4 stages of the selected
// coupling (PTF or CCF), the RK4 combination, and the coupling's post-step
// operation (PTF state relaxation, CCF correction taper).
void advanceL1Timestep(
    crescent::shared::LevelContext& level_context,
    amrex::Real time);

// Output derivative bundle of the current child state; no-op on levels
// other than the active child.
void recomputeL1OutputDerivatives(
    crescent::shared::LevelContext& level_context);

} // namespace crescent::nested::free_surface_evolution

#endif // CRESCENT_NESTED_TIMESTEP_H

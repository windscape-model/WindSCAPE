#ifndef CRESCENT_CORE_FREE_SURFACE_EVOLUTION_H
#define CRESCENT_CORE_FREE_SURFACE_EVOLUTION_H

#include <functional>

namespace crescent::shared {
class LevelContext;
}

namespace crescent::core::free_surface_evolution {

// post_stage_hook, when non-null, is invoked with the RK stage index after
// each stage's derivatives and RHS are final. It is observational (used by
// the nested correction-field coupling to capture parent stage data) and
// must not modify the core solution.
void advanceLevelTimestep(
    const crescent::shared::LevelContext& level_context,
    const std::function<void(int)>* post_stage_hook = nullptr);

} // namespace crescent::core::free_surface_evolution

#endif // CRESCENT_CORE_FREE_SURFACE_EVOLUTION_H

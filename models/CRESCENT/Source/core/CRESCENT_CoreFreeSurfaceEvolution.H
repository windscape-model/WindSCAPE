#ifndef CRESCENT_CORE_FREE_SURFACE_EVOLUTION_H
#define CRESCENT_CORE_FREE_SURFACE_EVOLUTION_H

namespace crescent::shared {
class LevelEvolutionContext;
}

namespace crescent::core::free_surface_evolution {

void advanceLevelTimestep(
    const crescent::shared::LevelEvolutionContext& level_context);

} // namespace crescent::core::free_surface_evolution

#endif // CRESCENT_CORE_FREE_SURFACE_EVOLUTION_H

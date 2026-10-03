// Primary nested-L1 timestep implementation.

#include "CRESCENT_NestedFreeSurfaceEvolution.H"
#include "CRESCENT_NestedRK4Stage.H"
#include "CRESCENT_NestedRuntimeLogging.H"
#include "CRESCENT_NestedStateBlend.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"

namespace crescent::nested::free_surface_evolution {

void advanceL1Timestep(crescent::shared::LevelEvolutionContext& level_context,
                       amrex::Real time)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    crescent::nested::runtime_logging::logNestedL1TimestepStart(
        solution.eta[lev],
        runtime.verbosity);

    const amrex::Real c[4] = {0.0, 0.5, 0.5, 1.0};

    for (int irk = 0; irk < 4; ++irk) {
        crescent::nested::rk4_stage::runL1RK4Stage(
            level_context,
            irk,
            time + c[irk] * runtime.dt[lev],
            c[irk]);
    }

    crescent::nested::state_blend::finalizeL1FreeSurfaceState(
        level_context);
}

} // namespace crescent::nested::free_surface_evolution

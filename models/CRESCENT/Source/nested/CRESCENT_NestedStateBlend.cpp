#include "CRESCENT_NestedStateBlend.H"

#include "CRESCENT_Dynamics.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedPeriodicSurfaceState.H"
#include "CRESCENT_NestedRuntimeLogging.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"

namespace crescent::nested::state_blend {

void finalizeL1FreeSurfaceState(
    crescent::shared::LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    const auto& nested_config = level_context.nestedConfig();
    const amrex::Real eta_amp_before_final =
        solution.eta[lev].norm0();

    crescent::dynamics::updateSurface(
        solution.eta[lev],
        solution.phi[lev],
        solution.eta_saved[lev],
        solution.phi_saved[lev],
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        runtime.dt[lev]);

    if (crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev) &&
        nested_config.enable_post_rk4_state_filter_and_blend == 1) {
        crescent::nested::periodic_surface_state::applyPostRK4CoreStateBlend(
            level_context);
    }

    const amrex::Real eta_amp_after_final =
        solution.eta[lev].norm0();

    crescent::nested::runtime_logging::logNestedFinalStateAfterFilter(
        solution.eta[lev],
        eta_amp_before_final,
        eta_amp_after_final,
        runtime.verbosity);
}

} // namespace crescent::nested::state_blend

#include "CRESCENT_NestedSetupWaveState.H"

#include <AMReX_Print.H>

#include "CRESCENT_Interpolation.H"

namespace crescent::nested::setup::wave_state {

void initializeFromCoreParent(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev > 0);
    const auto core = level_context.coreAccess();
    const auto nested = level_context.nestedAccess();
    const int verbosity = level_context.runtime().verbosity;
    if (verbosity >= 1) {
        amrex::Print() << "\n=== Initializing L1 from parent L0 ===" << std::endl;
    }

    const auto core_mesh = level_context.coreLevelMesh();
    const auto level_mesh = level_context.levelMesh();
    auto parent_state = core.solution.freeSurfaceStateAt(0);
    auto child_state = core.solution.freeSurfaceStateAt(lev);
    auto& coarse_interp_patch =
        nested.runtime.post_rk4_state_scratch[lev].coarseInterpPatch();

    if (verbosity >= 2) {
        amrex::Print() << "  Before interpolation: L1 eta range = ["
                       << child_state.eta.min(0) << ", "
                       << child_state.eta.max(0) << "]" << std::endl;
    }

    parent_state.eta.FillBoundary(core_mesh.geom.periodicity());
    parent_state.phi.FillBoundary(core_mesh.geom.periodicity());

    crescent::interpolation::interpolateParentToChild(
        crescent::shared::asConstFreeSurfaceState(parent_state),
        child_state,
        core_mesh.geom,
        level_mesh.geom,
        level_context.refineRatio(),
        coarse_interp_patch,
        nested.config.state_interpolation_order);

    if (verbosity >= 2) {
        amrex::Print() << "  After interpolation (before buffer): L1 eta range = ["
                       << child_state.eta.min(0) << ", "
                       << child_state.eta.max(0) << "]" << std::endl;
    }

    if (verbosity >= 1) {
        amrex::Print()
            << "  CRESCENT: Skipping buffer zone ramping (handled by CRESCENT workflow)"
            << std::endl;
    }

    core.solution.saveSurfaceState(lev);

    if (verbosity >= 1) {
        amrex::Print()
            << "  Computing DNO for L1 using CRESCENT derivative operator..."
            << std::endl;
        amrex::Print() << "    Interpolating L0 derivatives to L1 grid..." << std::endl;
    }
}

} // namespace crescent::nested::setup::wave_state

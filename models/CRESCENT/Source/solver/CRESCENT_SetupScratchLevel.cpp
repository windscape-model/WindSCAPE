#include "CRESCENT_SetupLifecycle.H"

#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

namespace {

void defineScratchLevelStateStorage(
    const crescent::solver::setup::ScratchLevelContext& scratch_level_context)
{
    const int lev = scratch_level_context.level();
    const auto core = scratch_level_context.coreAccess();
    const auto level_mesh = scratch_level_context.levelMesh();
    core.solution.defineScratchLevelStorage(
        lev,
        level_mesh.grids,
        level_mesh.dmap,
        core.wave.spectral.perturbation_order);
    core.root_layouts.defineLevel(lev, level_mesh.geom);
}

void initializeScratchLevelSurfaceState(
    const crescent::solver::setup::ScratchLevelContext& scratch_level_context)
{
    const int lev = scratch_level_context.level();
    const auto core = scratch_level_context.coreAccess();
    auto& problem = scratch_level_context.problem();
    const auto level_mesh = scratch_level_context.levelMesh();
    auto free_surface_state = core.solution.freeSurfaceStateAt(lev);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(free_surface_state.eta, amrex::TilingIfNotGPU());
         mfi.isValid();
         ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> eta_arr = free_surface_state.eta.array(mfi);
        amrex::Array4<amrex::Real> phi_arr = free_surface_state.phi.array(mfi);

        problem.init_data(
            bx, level_mesh.geom.data(), eta_arr, phi_arr);
    }
}

} // namespace

namespace crescent::solver::setup {

void initializeScratchLevelStorage(
    const ScratchLevelContext& scratch_level_context)
{
    defineScratchLevelStateStorage(scratch_level_context);
    initializeScratchLevelSurfaceState(scratch_level_context);
}

} // namespace crescent::solver::setup

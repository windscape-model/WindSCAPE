// ERF-CRESCENT coupling interface implemented on top of the CRESCENT solver.
//
// These methods are only exercised when CRESCENT is driven externally (e.g. by
// ERF). Standalone runs never call them, so standalone behaviour is unchanged.

#include "CRESCENT.H"

#include <AMReX_Print.H>

#include "CRESCENT_Checkpoint.H"
#include "CRESCENT_FreeSurfaceEvolution.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_SolverState.H"

void CRESCENT::setAtmosphericPressure(const amrex::MultiFab& p_atm, amrex::Real time)
{
    // Backward-compatible entry point: receive on the base level (the level coupled to ERF).
    setAtmosphericPressure(p_atm, time, 0);
}

void CRESCENT::setAtmosphericPressure(const amrex::MultiFab& p_atm, amrex::Real time, int lev)
{
    // Level-aware "receive" half of the external coupling exchange: overwrite the surface
    // pressure used in computeSurfaceRHS on level `lev` with the externally supplied (ERF)
    // forcing. `lev` lets a driver target CRESCENT's base level (current coupling) or, in the
    // future, a refined level in a region of interest.
    auto& pressure = solverState().core.solution.pressure;

    if (lev < 0 || lev >= static_cast<int>(pressure.size())) {
        amrex::Abort("CRESCENT::setAtmosphericPressure: requested level out of range");
    }
    if (pressure[lev].boxArray() != p_atm.boxArray()) {
        amrex::Abort("CRESCENT::setAtmosphericPressure: BoxArray mismatch between external pressure and CRESCENT surface field");
    }

    // Overwrite the surface pressure used in computeSurfaceRHS with the external forcing.
    amrex::MultiFab::Copy(pressure[lev], p_atm, 0, 0, 1, 0);
    use_atmospheric_pressure = true;

    const auto& runtime = solverState().core.runtime;
    if (runtime.verbosity >= 1) {
        amrex::Print() << "  [CRESCENT] Received atmospheric pressure at t=" << time
                       << " (level " << lev << ")"
                       << ", min=" << pressure[lev].min(0)
                       << ", max=" << pressure[lev].max(0) << std::endl;
    }
}

void CRESCENT::advanceSingleTimestep(amrex::Real /*dt*/)
{
    // Advance the free surface by one solver timestep (internal dt[0]), mirroring
    // run_loop::executeSingleStep minus the per-step reporting and file output.
    auto context = solverContext();
    amrex::Real time = solverState().core.runtime.t_new[0];
    crescent::solver::free_surface_evolution::advanceAllLevelsForStep(
        context, time);
    crescent::solver::free_surface_evolution::updateTimeStateAfterStep(
        context, time);
}

const amrex::MultiFab& CRESCENT::getEta(int lev) const
{
    return solverState().core.solution.eta[lev];
}

const amrex::MultiFab& CRESCENT::getU(int lev) const
{
    return solverState().core.solution.u[lev];
}

const amrex::MultiFab& CRESCENT::getV(int lev) const
{
    return solverState().core.solution.v[lev];
}

const amrex::MultiFab& CRESCENT::getVerticalVelocity(int lev) const
{
    return solverState().core.solution.phi_z[lev];
}

bool CRESCENT::shouldWriteCheckpoint() const
{
    const auto& runtime = solverState().core.runtime;
    return runtime.checkpoint.interval > 0 &&
           (runtime.istep[0] % runtime.checkpoint.interval == 0);
}

void CRESCENT::WriteCheckpointFile()
{
    crescent::solver::checkpoint::writeCheckpointFile(solverContext());
}

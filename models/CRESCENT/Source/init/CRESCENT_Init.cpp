#include "CRESCENT.H"
#include "solver/CRESCENT_Checkpoint.H"
#include "solver/CRESCENT_Setup.H"
#include "CRESCENT_SolverContext.H"

// Primary implementation file.

CRESCENT::~CRESCENT() = default;

void CRESCENT::Init() {
    auto& solver = solverState();

    if (!solver.core.runtime.checkpoint.restart_file.empty()) {
        crescent::solver::checkpoint::readCheckpointFile({
            solverState(),
            max_level,
            grids,
            dmap,
            geom,
            refRatio(),
            [this](
                int restart_finest_level,
                const amrex::Vector<amrex::BoxArray>& checkpoint_grids,
                const amrex::Vector<amrex::Real>& restart_time) {
                SetFinestLevel(restart_finest_level);
                for (int lev = 0; lev <= finest_level; ++lev) {
                    amrex::DistributionMapping level_dmap =
                        MakeDistributionMap(lev, checkpoint_grids[lev]);
                    MakeNewLevelFromScratch(
                        lev,
                        restart_time[lev],
                        checkpoint_grids[lev],
                        level_dmap);
                    SetBoxArray(lev, checkpoint_grids[lev]);
                    SetDistributionMap(lev, level_dmap);
                }
            }});
        return;
    }

    InitFromScratch(solver.core.runtime.start_time);
    crescent::solver::setup::initializeSolverState(solverContext());
}

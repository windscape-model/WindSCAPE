#include "CRESCENT.H"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_TagBox.H>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedDomains.H"
#include "CRESCENT_SetupLifecycle.H"
#include "CRESCENT_SolverOwnership.H"

void CRESCENT::ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time, int ngrow)
{
    (void)time;
    (void)ngrow;
    const auto solver = solverOwnership().solverAccess();

    if (lev == 0 &&
        crescent::nested::hasActiveNestedL1(
            solver.nested.config,
            max_level)) {
        const auto context =
            crescent::nested::domains::makeL1RefinementTaggingContext(
                solver.nested.config,
                geom[lev],
                solver.core.runtime.verbosity,
                refRatio(0),
                blockingFactor(1));
        crescent::nested::domains::tagL1RefinementRegion(tags, context);
    }
}

void CRESCENT::MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm)
{
    const auto solver = solverOwnership().solverAccess();
    if (solver.core.runtime.verbosity >= 1) {
        amrex::Print() << "MakeNewLevelFromCoarse called for level " << lev << std::endl;
    }
    MakeNewLevelFromScratch(lev, time, ba, dm);
}

void CRESCENT::RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm)
{
    (void)lev;
    (void)time;
    (void)ba;
    (void)dm;

    amrex::Abort(
        "CRESCENT does not support AMReX RemakeLevel/regridding after initialization; "
        "the spectral and nested workspaces require a fixed level layout.");
}

void CRESCENT::ClearLevel(int lev)
{
    (void)lev;

    amrex::Abort(
        "CRESCENT does not support AMReX ClearLevel after initialization; "
        "restart with a fixed hierarchy instead of dynamic regridding.");
}

void CRESCENT::MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm)
{
    (void)time;

    crescent::solver::setup::initializeScratchLevelStorage(
        crescent::solver::setup::makeScratchLevelContext(
            lev,
            solverOwnership().coreAccess(),
            solverOwnership().amrAccess(),
            geom[lev],
            ba,
            dm));
}

#include "CRESCENT_CoreFreeSurfaceEvolution.H"

#include <AMReX.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_CoreTimestepLogging.H"
#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_Dynamics.H"
#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_Spectral.H"

namespace crescent::core::free_surface_evolution {

namespace {

void runRKStage(
    const crescent::shared::LevelEvolutionContext& level_context,
    int rk_stage)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    auto& spectral_workspace = level_context.spectralWorkspace();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            level_context.levelMesh().domain());
    const auto& spectral_settings = level_wave_operator.spectral_settings;

    crescent::core::logging::logRKStageHeader(runtime, rk_stage);

    crescent::dynamics::advanceSurfaceRK4(
        solution.eta[lev],
        solution.phi[lev],
        solution.eta_saved[lev],
        solution.phi_saved[lev],
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        runtime.dt[lev],
        rk_stage);

    crescent::core::logging::logPostAdvanceSurfaceState(
        solution,
        runtime,
        lev,
        rk_stage);

    crescent::spectral::applyLowPassFilter(
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        solution.eta[lev],
        solution.phi[lev],
        level_wave_operator.spectral_workspace.eta_spectral[lev],
        spectral_settings,
        level_wave_operator.spectral_workspace.peak_wavenumber_gather_scratch[lev],
        level_wave_operator.spectral_workspace.peak_wavenumber_reduction_scratch[lev]);

    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        core,
        solution.derivativeBundleAt(lev));

    crescent::core::logging::logPreSurfaceRhsDiagnostics(
        solution,
        runtime,
        lev,
        rk_stage);

    if (lev != 0) {
        amrex::Abort("Core free-surface evolution supports only the pure L0 path");
    }

    crescent::dynamics::computeSurfaceRHS(
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        solution.eta[lev],
        solution.eta_x[lev],
        solution.eta_y[lev],
        solution.phi_x[lev],
        solution.phi_y[lev],
        solution.phi_z[lev],
        solution.pressure[lev],
        spectral_workspace.compute_rhs_metric[lev],
        spectral_workspace.compute_rhs_phi_z_sq[lev],
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        crescent::constants::INV_GRAVITY,
        level_wave_operator.perturbation_order,
        rk_stage,
        spectral_workspace.dealiasing_scratch[lev]);

    crescent::core::logging::logPostSurfaceRhsDiagnostics(
        solution,
        runtime,
        lev,
        rk_stage);
}

} // namespace

void advanceLevelTimestep(
    const crescent::shared::LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();

    runRKStage(level_context, 0);

    for (int rk_stage = 1; rk_stage < 4; ++rk_stage) {
        runRKStage(level_context, rk_stage);
    }

    crescent::dynamics::updateSurface(
        solution.eta[lev],
        solution.phi[lev],
        solution.eta_saved[lev],
        solution.phi_saved[lev],
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        runtime.dt[lev]);
}

} // namespace crescent::core::free_surface_evolution

#include "CRESCENT_CoreFreeSurfaceEvolution.H"

#include <AMReX.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_Dynamics.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_Spectral.H"

namespace crescent::core::free_surface_evolution {

namespace {

void runRKStage(
    const crescent::shared::LevelContext& level_context,
    int rk_stage)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    auto& spectral_workspace = level_context.spectralWorkspace();
    const auto level_wave_operator = level_context.waveOperator();
    const auto& spectral_settings = level_wave_operator.spectral_settings;

    crescent::dynamics::advanceSurfaceRK4(
        solution.eta[lev],
        solution.phi[lev],
        solution.eta_saved[lev],
        solution.phi_saved[lev],
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        runtime.dt[lev],
        rk_stage);

    crescent::spectral::applyLowPassFilter(
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        solution.eta[lev],
        solution.phi[lev],
        spectral_settings,
        spectral_workspace.filter_weights[lev]);

    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        core,
        solution.derivativeBundleAt(lev));

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

    if (core.wave.spectral.adjustment_periods > 0.0) {
        // Dommermuth adjustment: t_new[lev] holds the step-start time during
        // the advance; the ramp is applied before the nested capture hooks so
        // interpolated core tendencies stay ramp-consistent on the child.
        constexpr amrex::Real stage_c[4] = {0.0, 0.5, 0.5, 1.0};
        const amrex::Real t_peak = 2.0 * crescent::constants::PI /
            core.wave.initialization.peak_frequency;
        crescent::dynamics::applyAdjustmentRampToStageRHS(
            solution.rhs_eta[lev],
            solution.rhs_phi[lev],
            solution.eta[lev],
            solution.phi[lev],
            solution.pressure[lev],
            level_wave_operator.fft,
            level_wave_operator.geomdata,
            runtime.t_new[lev] + stage_c[rk_stage] * runtime.dt[lev],
            core.wave.spectral.adjustment_periods * t_peak,
            rk_stage,
            spectral_workspace.dealiasing_scratch[lev]);
    }
}

} // namespace

void advanceLevelTimestep(
    const crescent::shared::LevelContext& level_context,
    const std::function<void(int)>* post_stage_hook)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();

    runRKStage(level_context, 0);
    if (post_stage_hook != nullptr && *post_stage_hook) {
        (*post_stage_hook)(0);
    }

    for (int rk_stage = 1; rk_stage < 4; ++rk_stage) {
        runRKStage(level_context, rk_stage);
        if (post_stage_hook != nullptr && *post_stage_hook) {
            (*post_stage_hook)(rk_stage);
        }
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

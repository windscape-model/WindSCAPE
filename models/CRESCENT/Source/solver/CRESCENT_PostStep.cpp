#include "CRESCENT_PostStep.H"

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_Spectral.H"

namespace crescent::solver::poststep {

void applyOptionalPostTimestepLowPassFilter(
    const crescent::shared::LevelEvolutionContext& level_context)
{
    if (level_context.nestedConfig().enable_post_timestep_low_pass_filter != 1) {
        return;
    }

    const int lev = level_context.level();
    auto free_surface_state = level_context.solutionState().freeSurfaceStateAt(lev);
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());
    const auto& spectral_settings = level_wave_operator.spectral_settings;

    crescent::spectral::applyLowPassFilter(
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        free_surface_state,
        level_wave_operator.spectral_workspace.eta_spectral[lev],
        spectral_settings,
        level_wave_operator.spectral_workspace.peak_wavenumber_gather_scratch[lev],
        level_wave_operator.spectral_workspace.peak_wavenumber_reduction_scratch[lev]);
}

void recomputeStandardDerivativeBundle(
    const crescent::shared::LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());
    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        level_context.coreAccess(),
        level_context.solutionState().derivativeBundleAt(lev));
}

void computeVelocityOutputs(
    const crescent::shared::LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());
    auto derivative_bundle =
        crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev)
            ? level_context.nestedRuntime().coupling_derivatives.l1_blended_derivatives.bundleAt(lev)
            : solution.derivativeBundleAt(lev);

    crescent::physics::computeHorizontalVelocities(
        solution.u[lev],
        solution.v[lev],
        derivative_bundle,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev]);
}

} // namespace crescent::solver::poststep

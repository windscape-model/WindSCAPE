#include "CRESCENT_PostStep.H"

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_Spectral.H"

namespace crescent::solver::poststep {

void recomputeStandardDerivativeBundle(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    const auto level_wave_operator = level_context.waveOperator();
    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        level_context.core(),
        level_context.solutionState().derivativeBundleAt(lev));
}

void computeVelocityOutputs(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto level_wave_operator = level_context.waveOperator();
    const bool is_nested_l1 =
        crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(), lev);
    auto derivative_bundle = is_nested_l1
            ? level_context.nestedRuntime().coupling_derivatives.l1_blended_derivatives.bundleAt(lev)
            : solution.derivativeBundleAt(lev);
    // The blended child bundle is non-periodic: no spectral dealiasing.
    const bool apply_spectral_dealiasing = !is_nested_l1;

    crescent::physics::computeHorizontalVelocities(
        solution.u[lev],
        solution.v[lev],
        derivative_bundle,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev],
        apply_spectral_dealiasing);
}

} // namespace crescent::solver::poststep

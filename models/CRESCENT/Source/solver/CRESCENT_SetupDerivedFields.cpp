#include "CRESCENT_SetupDerivedFields.H"

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_PhysicsInitializationLogging.H"
#include "CRESCENT_Physics.H"

namespace crescent::solver::setup::derived_fields {

void initializeCoreDerivedFields(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev == 0);
    const auto core = level_context.coreAccess();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            level_context.levelMesh().domain());
    auto derivative_bundle = core.solution.derivativeBundleAt(lev);
    auto& solution = core.solution;

    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        core,
        derivative_bundle);

    crescent::physics::initialization_logging::logSurfaceDerivativeSummary(
        solution.eta[lev],
        solution.eta_x[lev],
        solution.eta_y[lev],
        core.runtime.verbosity);
    crescent::physics::initialization_logging::logPotentialSpectrum(
        solution.phi[lev],
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        core.runtime.verbosity);
    crescent::physics::initialization_logging::logPhiExpansionState(
        solution.phi_expansion[lev],
        solution.eta_pow[lev],
        level_wave_operator.perturbation_order,
        level_wave_operator.geomdata,
        core.runtime.verbosity);
    crescent::physics::initialization_logging::logVerticalVelocityState(
        solution.phi_z[lev],
        level_wave_operator.geomdata,
        core.runtime.verbosity);

    crescent::physics::computeHorizontalVelocities(
        solution.u[lev],
        solution.v[lev],
        derivative_bundle,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev]);
}

} // namespace crescent::solver::setup::derived_fields

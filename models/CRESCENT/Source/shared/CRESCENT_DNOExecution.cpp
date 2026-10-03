#include "CRESCENT_DNOExecution.H"

#include "CRESCENT_DNOConfig.H"
#include "CRESCENT_DNOWorkspaceViews.H"
#include "CRESCENT_NestedDNO.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"

namespace crescent::dno {

SurfaceStateDNOContext makeSurfaceStateDNOContext(
    crescent::shared::LevelWaveOperatorView level_wave_operator,
    int dno_verbosity,
    crescent::shared::ConstFreeSurfaceStateView source_state,
    const amrex::MultiFab& wavenumber_pow_source,
    amrex::MultiFab& eta_pow_storage,
    amrex::MultiFab& phi_expansion_storage,
    amrex::MultiFab& phi_z_out)
{
    return {
        level_wave_operator,
        dno_verbosity,
        source_state,
        wavenumber_pow_source,
        eta_pow_storage,
        phi_expansion_storage,
        phi_z_out
    };
}

} // namespace crescent::dno

namespace {

struct StandardDerivativeBundleContext
{
    crescent::dno::SurfaceStateDNOContext surface_state;
    crescent::shared::DerivativeBundleView derivative_bundle;
};

struct PreparedSurfaceStateDNOExecution
{
    crescent::dno::WorkspaceView workspace;
    crescent::config::DNOConfig config;
};

crescent::dno::SurfaceStateDNOContext makeSolutionStateDNOContext(
    crescent::shared::LevelWaveOperatorView level_wave_operator,
    const crescent::ownership::CoreSolverAccess& core,
    amrex::MultiFab& phi_z_out)
{
    const int lev = level_wave_operator.lev;
    return crescent::dno::makeSurfaceStateDNOContext(
        level_wave_operator,
        core.runtime.verbosity,
        crescent::shared::asConstFreeSurfaceState(
            core.solution.freeSurfaceStateAt(lev)),
        core.solution.wavenumber_pow[lev],
        core.solution.eta_pow[lev],
        core.solution.phi_expansion[lev],
        phi_z_out);
}

StandardDerivativeBundleContext makeSolutionStateStandardDerivativeBundleContext(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const crescent::ownership::CoreSolverAccess& core,
    crescent::shared::DerivativeBundleView derivative_bundle)
{
    return {
        makeSolutionStateDNOContext(
            level_wave_operator,
            core,
            derivative_bundle.phi_z),
        derivative_bundle
    };
}

crescent::dno::WorkspaceView makeDNOWorkspace(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    amrex::MultiFab& eta_pow_storage,
    amrex::MultiFab& phi_expansion_storage)
{
    const int lev = level_wave_operator.lev;
    auto& spectral_workspace = level_wave_operator.spectral_workspace;
    return crescent::dno::WorkspaceView::fromExisting(
        eta_pow_storage, phi_expansion_storage,
        spectral_workspace.expansion_deriv[lev], spectral_workspace.deriv_spectral[lev],
        spectral_workspace.expansion_accum[lev], spectral_workspace.spectral_work[lev],
        spectral_workspace.phi_expansion_spectrum[lev]);
}

PreparedSurfaceStateDNOExecution prepareSurfaceStateDNOExecution(
    const crescent::dno::SurfaceStateDNOContext& context)
{
    return {
        makeDNOWorkspace(
            context.level_wave_operator,
            context.eta_pow_storage,
            context.phi_expansion_storage),
        crescent::config::makeDNOConfig(
            context.level_wave_operator.perturbation_order,
            context.dno_verbosity)
    };
}

template <typename ExpandVelocityPotentialFn, typename ComputeVerticalVelocityFn>
void computeSurfaceStateDNO(
    const crescent::dno::SurfaceStateDNOContext& context,
    ExpandVelocityPotentialFn expand_velocity_potential,
    ComputeVerticalVelocityFn compute_vertical_velocity)
{
    const auto& level_wave_operator = context.level_wave_operator;
    auto& dealiasing_scratch =
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev];

    crescent::physics::computeEtaPowers(
        context.eta_pow_storage,
        context.source_state.eta,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.perturbation_order,
        dealiasing_scratch);

    auto prepared_execution = prepareSurfaceStateDNOExecution(context);
    const auto phi_expansion_workspace =
        crescent::dno::makePhiExpansionWorkspaceView(prepared_execution.workspace);
    const auto vertical_velocity_workspace =
        crescent::dno::makeVerticalVelocityWorkspaceView(prepared_execution.workspace);

    expand_velocity_potential(
        context.source_state.phi,
        context.wavenumber_pow_source,
        phi_expansion_workspace,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        prepared_execution.config);

    compute_vertical_velocity(
        context.phi_z_out,
        context.wavenumber_pow_source,
        vertical_velocity_workspace,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        prepared_execution.config,
        dealiasing_scratch);
}

void computeStandardDerivativeDNOBundle(
    const StandardDerivativeBundleContext& context)
{
    const auto& surface_state = context.surface_state;
    const auto& derivative_bundle = context.derivative_bundle;
    const auto& level_wave_operator = surface_state.level_wave_operator;
    const auto& wavenumbers = level_wave_operator.wavenumbers();

    crescent::physics::computeSurfaceDerivatives(
        derivative_bundle.eta_x,
        derivative_bundle.eta_y,
        derivative_bundle.phi_x,
        derivative_bundle.phi_y,
        surface_state.source_state.eta,
        surface_state.source_state.phi,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        level_wave_operator.spectral_workspace.directional_input_spectral[level_wave_operator.lev],
        level_wave_operator.spectral_workspace.directional_x_spectral[level_wave_operator.lev],
        level_wave_operator.spectral_workspace.directional_y_spectral[level_wave_operator.lev]);

    computeSurfaceStateDNO(
        surface_state,
        crescent::physics::expandVelocityPotential,
        crescent::physics::computeVerticalVelocity);
}

} // namespace

namespace crescent::dno {

void computeSolutionStateDNO(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const crescent::ownership::CoreSolverAccess& core,
    amrex::MultiFab& phi_z_out)
{
    computeSurfaceStateDNO(
        ::makeSolutionStateDNOContext(
            level_wave_operator,
            core,
            phi_z_out),
        crescent::physics::expandVelocityPotential,
        crescent::physics::computeVerticalVelocity);
}

void computeSolutionStateDerivativeBundle(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const crescent::ownership::CoreSolverAccess& core,
    crescent::shared::DerivativeBundleView derivative_bundle)
{
    computeStandardDerivativeDNOBundle(
        makeSolutionStateStandardDerivativeBundleContext(
            level_wave_operator,
            core,
            derivative_bundle));
}

void computeNestedL1DNO(const SurfaceStateDNOContext& context)
{
    computeSurfaceStateDNO(
        context,
        crescent::nested::derivatives::expandVelocityPotential,
        crescent::nested::derivatives::computeVerticalVelocity);
}

} // namespace crescent::dno

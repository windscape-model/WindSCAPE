#include "CRESCENT_DNOExecution.H"

#include "CRESCENT_Constants.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_SolverState.H"
#include "CRESCENT_SpectralGridMetrics.H"

namespace crescent::dno {

namespace {

crescent::physics::VerticalVelocityWorkspace levelWorkspace(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    amrex::MultiFab& eta_pow,
    amrex::MultiFab& phi_expansion)
{
    const int lev = level_wave_operator.lev;
    auto& workspace = level_wave_operator.spectral_workspace;
    return {
        eta_pow,
        phi_expansion,
        workspace.expansion_deriv[lev],
        workspace.deriv_spectral[lev],
        workspace.expansion_accum[lev],
        workspace.spectral_work[lev],
        workspace.phi_expansion_spectrum[lev],
        workspace.dealiasing_scratch[lev],
        workspace.vertical_velocity_kpow[lev]
    };
}

} // namespace

void computeSolutionStateDerivativeBundle(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    crescent::CoreSolverState& core,
    crescent::shared::DerivativeBundleView derivative_bundle)
{
    const int lev = level_wave_operator.lev;
    auto& solution = core.solution;
    auto& workspace = level_wave_operator.spectral_workspace;
    const auto& wavenumbers = level_wave_operator.wavenumbers();
    const auto& geom = level_wave_operator.geomdata;

    crescent::physics::computeSurfaceDerivatives(
        derivative_bundle.eta_x,
        derivative_bundle.eta_y,
        derivative_bundle.phi_x,
        derivative_bundle.phi_y,
        solution.eta[lev],
        solution.phi[lev],
        level_wave_operator.fft,
        geom,
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        workspace.directional_input_spectral[lev],
        workspace.directional_x_spectral[lev],
        workspace.directional_y_spectral[lev]);

    const amrex::Real Lx = geom.ProbHi()[0] - geom.ProbLo()[0];
    const amrex::Real Ly = geom.ProbHi()[1] - geom.ProbLo()[1];
    crescent::physics::computeVerticalVelocity(
        derivative_bundle.phi_z,
        solution.eta[lev],
        solution.phi[lev],
        solution.wavenumber_pow[lev],
        levelWorkspace(
            level_wave_operator,
            solution.eta_pow[lev],
            solution.phi_expansion[lev]),
        level_wave_operator.fft,
        geom,
        level_wave_operator.perturbation_order,
        crescent::constants::PI2 / Lx,
        crescent::constants::PI2 / Ly,
        true);
}

void computeNestedL1DNO(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const amrex::MultiFab& eta,
    const amrex::MultiFab& phi,
    const amrex::MultiFab& wavenumber_pow,
    amrex::MultiFab& eta_pow,
    amrex::MultiFab& phi_expansion,
    amrex::MultiFab& phi_z,
    bool apply_final_projection)
{
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(
            phi_expansion,
            level_wave_operator.geomdata);
    crescent::physics::computeVerticalVelocity(
        phi_z,
        eta,
        phi,
        wavenumber_pow,
        levelWorkspace(level_wave_operator, eta_pow, phi_expansion),
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.perturbation_order,
        spectral_grid.dk_x,
        spectral_grid.dk_y,
        apply_final_projection);
}

} // namespace crescent::dno

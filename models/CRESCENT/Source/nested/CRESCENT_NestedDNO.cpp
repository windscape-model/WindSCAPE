#include "CRESCENT_NestedDNO.H"

#include "CRESCENT_NestedDerivativeDiagnostics.H"
#include "CRESCENT_DNO_Internal.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_PerformanceTimer.H"

namespace crescent::nested::derivatives {

using amrex::cMultiFab;
using amrex::GeometryData;
using amrex::MultiFab;
using amrex::Real;

namespace FFT = amrex::FFT;

// The earlier Ω_S-only derivative path was retired because it used the wrong
// spectral spacing (2π/Lx_S). The live nested solver computes spectral work on
// the full L1 grid only, with dk = 2π/Lx_L1.

void expandVelocityPotential(
    const MultiFab& phi,
    const MultiFab& wavenumber_pow,
    crescent::dno::PhiExpansionWorkspaceView workspace,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom,
    const crescent::config::DNOConfig& config)
{
    const int perturbation_order = config.perturbation_order;

    const auto execution_hooks =
        crescent::dno::detail::makePhiExpansionExecutionHooks(
            [](const cMultiFab&) {},
            [](const cMultiFab&, int, int) {},
            [](const MultiFab&, const MultiFab&, const GeometryData&, int, int) {});

    crescent::dno::detail::executePhiExpansionRecurrence(
        phi,
        wavenumber_pow,
        workspace,
        perturbation_order,
        fft,
        geom,
        execution_hooks);
}

void computeVerticalVelocity(
    MultiFab& phi_z,
    const MultiFab& wavenumber_pow,
    crescent::dno::VerticalVelocityWorkspaceView workspace,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom,
    const crescent::config::DNOConfig& config,
    MultiFab& dealiasing_scratch)
{
    // L1 path does not dealias inside the DNO kernel; the dealiasing scratch
    // argument is accepted for signature parity with the L0 path so the two
    // implementations can share the computeSurfaceStateDNO dispatch template.
    (void)dealiasing_scratch;

    crescent::diagnostics::ScopedTimer timer(
        "computeVerticalVelocity_L1", config.verbosity);

    const int perturbation_order = config.perturbation_order;
    const int verbosity = config.verbosity;
    const auto diagnostics =
        crescent::nested::diagnostics::makeVerticalVelocityDiagnosticHooks(
            perturbation_order,
            verbosity);
    const auto execution_hooks = diagnostics.makeExecutionHooks();

    diagnostics.logInputs(
        workspace.phi_expansion,
        wavenumber_pow);

    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(
            workspace.phi_expansion,
            geom);
    crescent::dno::detail::executeVerticalVelocityAccumulation(
        phi_z,
        workspace,
        perturbation_order,
        fft,
        spectral_grid.dk_x,
        spectral_grid.dk_y,
        execution_hooks);

    diagnostics.logSuccess(phi_z);
}

} // namespace crescent::nested::derivatives

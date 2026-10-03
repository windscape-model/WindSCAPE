#include "CRESCENT_Physics.H"
#include "CRESCENT_DNO_Internal.H"
#include "CRESCENT_DNOWorkspaceViews.H"
#include "CRESCENT_PhysicsExecutionDiagnostics.H"
#include "CRESCENT_Physics_Helpers.H"
#include "CRESCENT_SpectralDerivativeExecution.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_PerformanceTimer.H"
#include <AMReX_GpuComplex.H>

namespace {

using amrex::cMultiFab;
using amrex::GeometryData;
using amrex::MultiFab;
using amrex::Real;

namespace FFT = amrex::FFT;

struct DomainWavenumberSpacing
{
    Real dk_x;
    Real dk_y;
};

DomainWavenumberSpacing computeDomainWavenumberSpacing(const GeometryData& geom)
{
    Real Lx = geom.ProbHi()[0] - geom.ProbLo()[0];
    Real Ly = geom.ProbHi()[1] - geom.ProbLo()[1];
    return {crescent::constants::PI2 / Lx, crescent::constants::PI2 / Ly};
}

void computeDirectionalDerivativesNoFilter(
    MultiFab& field_x,
    MultiFab& field_y,
    const MultiFab& field,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    Real dkx,
    Real dky,
    cMultiFab& input_spectral_scratch,
    cMultiFab& x_spectral_scratch,
    cMultiFab& y_spectral_scratch)
{
    const auto unit_spectral_weight =
        [] AMREX_GPU_HOST_DEVICE(int, int) noexcept -> Real {
            return Real(1.0);
        };

    crescent::spectral::executeDirectionalDerivativePair(
        field_x,
        field_y,
        field,
        fft,
        spectral_grid,
        dkx,
        dky,
        unit_spectral_weight,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
}

} // namespace

namespace crescent::physics {

using amrex::Box;
using amrex::cMultiFab;
using amrex::GeometryData;
using amrex::make_alias;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;

namespace FFT = amrex::FFT;

void computeEtaPowers(MultiFab& eta_pow,
                      const MultiFab& eta,
                      std::unique_ptr<FFT::R2C<Real>>& fft,
                      const GeometryData& geom,
                      const int perturbation_order,
                      MultiFab& dealiasing_scratch)
{
    // Use helper function for eta power computation
    helpers::computeEtaPowersWithFactorial(eta_pow, eta, perturbation_order);

    // Dealias each component after computation
    for (int k = 1; k < perturbation_order - 1; ++k)
    {
        MultiFab tmp(eta_pow, make_alias, k, 1);
        crescent::spectral::dealiasing(tmp, fft, geom, dealiasing_scratch);
    }
}

void computeSurfaceDerivatives(MultiFab& eta_x, MultiFab& eta_y,
                                MultiFab& phi_x, MultiFab& phi_y,
                                const MultiFab& eta,
                                const MultiFab& phi,
                                std::unique_ptr<FFT::R2C<Real>>& fft,
                                const GeometryData& geom,
                                const Real dk_x, const Real dk_y,
                                cMultiFab& input_spectral_scratch,
                                cMultiFab& x_spectral_scratch,
                                cMultiFab& y_spectral_scratch)
{
    // Compute each field's x/y derivatives from one forward FFT, instead of
    // separate x and y transforms for the same input field.

    const auto spectral_grid = crescent::shared::makeSpectralGridMetrics(eta, geom);

    computeDirectionalDerivativesNoFilter(
        eta_x,
        eta_y,
        eta,
        fft,
        spectral_grid,
        dk_x,
        dk_y,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
    computeDirectionalDerivativesNoFilter(
        phi_x,
        phi_y,
        phi,
        fft,
        spectral_grid,
        dk_x,
        dk_y,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
}

void expandVelocityPotential(const MultiFab& phi,
                              const MultiFab& wavenumber_pow,
                              crescent::dno::PhiExpansionWorkspaceView workspace,
                              std::unique_ptr<FFT::R2C<Real>>& fft,
                              const GeometryData& geom,
                              const crescent::config::DNOConfig& config)
{
    crescent::diagnostics::ScopedTimer timer(
        "expandVelocityPotential", config.verbosity);

    const int perturbation_order = config.perturbation_order;
    const int verbosity = config.verbosity;
    const auto diagnostics =
        crescent::physics::execution_diagnostics::makePhiExpansionDiagnosticHooks(
            verbosity);

    crescent::dno::detail::executePhiExpansionRecurrence(
        phi,
        wavenumber_pow,
        workspace,
        perturbation_order,
        fft,
        geom,
        diagnostics.makeExecutionHooks());
}

void computeVerticalVelocity(MultiFab& phi_z,
                              const MultiFab& wavenumber_pow,
                              crescent::dno::VerticalVelocityWorkspaceView workspace,
                              std::unique_ptr<FFT::R2C<Real>>& fft,
                              const GeometryData& geom,
                              const crescent::config::DNOConfig& config,
                              MultiFab& dealiasing_scratch)
{
    crescent::diagnostics::ScopedTimer timer(
        "computeVerticalVelocity", config.verbosity);

    const int perturbation_order = config.perturbation_order;
    const int verbosity = config.verbosity;
    const auto diagnostics =
        crescent::physics::execution_diagnostics::makeVerticalVelocityDiagnosticHooks(
            perturbation_order,
            verbosity);
    const auto execution_hooks = diagnostics.makeExecutionHooks();

    diagnostics.logInputs(
        workspace.phi_expansion,
        workspace.eta_pow,
        wavenumber_pow);

    const DomainWavenumberSpacing domain_wavenumbers =
        computeDomainWavenumberSpacing(geom);
    crescent::dno::detail::executeVerticalVelocityAccumulation(
        phi_z,
        workspace,
        perturbation_order,
        fft,
        domain_wavenumbers.dk_x,
        domain_wavenumbers.dk_y,
        execution_hooks);

    diagnostics.logPreDealiasing(phi_z);

    crescent::spectral::dealiasing(phi_z, fft, geom, dealiasing_scratch);

    diagnostics.logPostDealiasing(phi_z);
}

void computeHorizontalVelocities(MultiFab& u, MultiFab& v,
                                  const MultiFab& phi_z,
                                  const MultiFab& eta_x, const MultiFab& eta_y,
                                  const MultiFab& phi_x, const MultiFab& phi_y,
                                  std::unique_ptr<FFT::R2C<Real>>& fft,
                                  const GeometryData& geom,
                                  MultiFab& dealiasing_scratch)
{
    for (MFIter mfi(eta_x, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto us_arr = u.array(mfi);
        auto vs_arr = v.array(mfi);
        auto phi_z_arr = phi_z.const_array(mfi);
        auto eta_x_arr = eta_x.const_array(mfi);
        auto eta_y_arr = eta_y.const_array(mfi);
        auto phi_x_arr = phi_x.const_array(mfi);
        auto phi_y_arr = phi_y.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // For collocated grids (non-staggered), directly compute the surface velocities
            us_arr(i,j,k) = phi_x_arr(i,j,k) - phi_z_arr(i,j,k) * eta_x_arr(i,j,k);
            vs_arr(i,j,k) = phi_y_arr(i,j,k) - phi_z_arr(i,j,k) * eta_y_arr(i,j,k);

            // Note: If moving to staggered grids in the future, interpolation would be needed here
            // For example, if u is staggered in x: us_arr(i,j,k) would need values from i-1/2 and i+1/2
        });
    }

    crescent::spectral::dealiasing(u, fft, geom, dealiasing_scratch);
    crescent::spectral::dealiasing(v, fft, geom, dealiasing_scratch);
}

void computeHorizontalVelocities(
    MultiFab& u,
    MultiFab& v,
    crescent::shared::DerivativeBundleView derivatives,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geom,
    MultiFab& dealiasing_scratch)
{
    computeHorizontalVelocities(
        u,
        v,
        derivatives.phi_z,
        derivatives.eta_x,
        derivatives.eta_y,
        derivatives.phi_x,
        derivatives.phi_y,
        fft,
        geom,
        dealiasing_scratch);
}

} // namespace crescent::physics

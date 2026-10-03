#include "CRESCENT_NestedSurfaceDerivatives.H"

#include "CRESCENT_NestedDerivativeDiagnostics.H"
#include "CRESCENT_SpectralDerivativeExecution.H"
#include "CRESCENT_SpectralFilterKernels.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_Constants.H"

#include <cmath>

using namespace amrex;

namespace crescent::nested::derivatives {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real linearDerivativeSpectralSigma(
    int filter_mode,
    int i,
    int ky_index,
    Real dk_x,
    Real dk_y,
    Real k_nyquist,
    Real k_peak,
    Real beta1,
    Real beta2);

void computeDirectionalDerivativesWithL1Filter(
    MultiFab& field_x_out,
    MultiFab& field_y_out,
    const MultiFab& field,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    int spectral_filter_mode,
    Real k_peak,
    Real beta1,
    Real beta2,
    cMultiFab& input_spectral_scratch,
    cMultiFab& x_spectral_scratch,
    cMultiFab& y_spectral_scratch,
    MultiFab& x_valid_scratch,
    MultiFab& y_valid_scratch)
{
    using namespace crescent::constants;
    const Real k_nyquist = spectral_grid.minNyquist();

    const auto spectral_weight =
        [=] AMREX_GPU_HOST_DEVICE(int i, int ky_index) noexcept -> Real {
            const int kx =
                crescent::spectral::filter_kernels::spectralKxIndex(
                    i,
                    spectral_grid.nx);
            const int cutoff_x =
                static_cast<int>(spectral_grid.nx * ONE_THIRD_CUTOFF);
            const int cutoff_y =
                static_cast<int>(spectral_grid.ny * ONE_THIRD_CUTOFF);

            if (crescent::spectral::filter_kernels::spectralModeExceedsCutoff(
                    kx,
                    ky_index,
                    cutoff_x,
                    cutoff_y)) {
                return Real(0.0);
            }

            return linearDerivativeSpectralSigma(
                spectral_filter_mode,
                i,
                ky_index,
                spectral_grid.dk_x,
                spectral_grid.dk_y,
                k_nyquist,
                k_peak,
                beta1,
                beta2);
        };

    crescent::spectral::executeDirectionalDerivativePairPreservingGrowCells(
        field_x_out,
        field_y_out,
        field,
        fft,
        spectral_grid,
        spectral_grid.dk_x,
        spectral_grid.dk_y,
        spectral_weight,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch,
        x_valid_scratch,
        y_valid_scratch);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real linearDerivativeSpectralSigma(
    int filter_mode,
    int i,
    int ky_index,
    Real dk_x,
    Real dk_y,
    Real k_nyquist,
    Real k_peak,
    Real beta1,
    Real beta2)
{
    if (filter_mode == 0) {
        return 1.0;
    }

    const Real k_mag =
        crescent::spectral::filter_kernels::spectralMagnitudeFromIndices(
            i,
            ky_index,
            dk_x,
            dk_y);
    Real sigma_total =
        crescent::spectral::filter_kernels::houLiWeight(
            k_mag,
            k_nyquist,
            crescent::constants::HOULIGHER_FILTER_ALPHA_AGGRESSIVE,
            static_cast<int>(crescent::constants::HOULIGHER_FILTER_POWER));

    if (filter_mode == 2 && k_peak > 0.0 && k_mag > k_peak) {
        sigma_total *= crescent::spectral::filter_kernels::lowPassWeight(
            k_mag,
            k_peak,
            beta1,
            beta2,
            true);
    }

    return sigma_total;
}

} // namespace

void computeSurfaceDerivatives(
    MultiFab& eta_x,
    MultiFab& eta_y,
    MultiFab& phi_x,
    MultiFab& phi_y,
    const MultiFab& eta,
    const MultiFab& phi,
    const Box& box_S,
    const Geometry& geom,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    int spectral_filter_mode,
    Real k_peak,
    Real beta1,
    Real beta2,
    cMultiFab& input_spectral_scratch,
    cMultiFab& x_spectral_scratch,
    cMultiFab& y_spectral_scratch,
    MultiFab& x_valid_scratch,
    MultiFab& y_valid_scratch,
    int verbosity)
{
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(
            eta,
            geom);

    crescent::nested::diagnostics::logSurfaceDerivativeComputationStart(
        spectral_grid.nx,
        spectral_grid.ny,
        spectral_grid.Lx,
        spectral_grid.Ly,
        spectral_grid.dk_x,
        spectral_grid.dk_y,
        spectral_filter_mode,
        box_S,
        verbosity);

    computeDirectionalDerivativesWithL1Filter(
        eta_x, eta_y, eta, fft,
        spectral_grid,
        spectral_filter_mode, k_peak, beta1, beta2,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch,
        x_valid_scratch,
        y_valid_scratch);
    crescent::nested::diagnostics::logSurfaceDerivativeRanges(
        "eta", eta_x, eta_y, verbosity);

    computeDirectionalDerivativesWithL1Filter(
        phi_x, phi_y, phi, fft,
        spectral_grid,
        spectral_filter_mode, k_peak, beta1, beta2,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch,
        x_valid_scratch,
        y_valid_scratch);
    crescent::nested::diagnostics::logSurfaceDerivativeRanges(
        "phi", phi_x, phi_y, verbosity);

    crescent::nested::diagnostics::logSurfaceDerivativeComputationEnd(verbosity);
}

} // namespace crescent::nested::derivatives

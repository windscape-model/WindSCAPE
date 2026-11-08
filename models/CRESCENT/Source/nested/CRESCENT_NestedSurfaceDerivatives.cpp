#include "CRESCENT_NestedSurfaceDerivatives.H"

#include "CRESCENT_Constants.H"
#include "CRESCENT_SpectralDerivativeExecution.H"
#include "CRESCENT_SpectralFilterKernels.H"
#include "CRESCENT_SpectralGridMetrics.H"

using namespace amrex;

namespace crescent::nested::derivatives {

namespace {

void computeDirectionalDerivatives(
    MultiFab& field_x_out,
    MultiFab& field_y_out,
    const MultiFab& field,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    cMultiFab& input_spectral_scratch,
    cMultiFab& x_spectral_scratch,
    cMultiFab& y_spectral_scratch,
    MultiFab& x_valid_scratch,
    MultiFab& y_valid_scratch)
{
    using namespace crescent::constants;

    // Weight 1 inside the square two-thirds cut, 0 outside.
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
            return Real(1.0);
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

} // namespace

void computeSurfaceDerivatives(
    MultiFab& eta_x,
    MultiFab& eta_y,
    MultiFab& phi_x,
    MultiFab& phi_y,
    const MultiFab& eta,
    const MultiFab& phi,
    const Geometry& geom,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    cMultiFab& input_spectral_scratch,
    cMultiFab& x_spectral_scratch,
    cMultiFab& y_spectral_scratch,
    MultiFab& x_valid_scratch,
    MultiFab& y_valid_scratch)
{
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(eta, geom);

    computeDirectionalDerivatives(
        eta_x, eta_y, eta, fft, spectral_grid,
        input_spectral_scratch, x_spectral_scratch, y_spectral_scratch,
        x_valid_scratch, y_valid_scratch);
    computeDirectionalDerivatives(
        phi_x, phi_y, phi, fft, spectral_grid,
        input_spectral_scratch, x_spectral_scratch, y_spectral_scratch,
        x_valid_scratch, y_valid_scratch);
}

} // namespace crescent::nested::derivatives

#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_SpectralFilterExecution.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_SpectralFilterKernels.H"

// Active nested-L1 filter implementation.

namespace crescent::nested::filter {

namespace FFT = amrex::FFT;

using amrex::GeometryData;
using amrex::GpuComplex;
using amrex::MultiFab;
using amrex::Real;

void applyHouLiFilter_L1(
    MultiFab& eta,
    MultiFab& phi,
    std::unique_ptr<FFT::R2C<Real>>& fft,
    const GeometryData& geomdata,
    const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
    crescent::shared::SpectralWeightTable& weights)
{
    const Real dk_x = spectral_settings.wavenumbers.dk_x;
    const Real dk_y = spectral_settings.wavenumbers.dk_y;
    const Real k_peak = spectral_settings.filter.peak_wavenumber_physical;
    const Real beta1 = spectral_settings.filter.beta1_filter;
    const Real beta2 = spectral_settings.filter.beta2_filter;
    const Real rolloff_alpha = spectral_settings.filter.l1_rolloff_alpha;
    const int rolloff_power = spectral_settings.filter.l1_rolloff_power;

    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(eta, geomdata);
    const int nx = spectral_grid.nx;
    const int ny = spectral_grid.ny;

    // k_nyquist from the child grid (nx/2 modes at spacing dk_x); the level
    // geometry describes the parent domain.
    const Real k_nyquist = spectral_grid.k_nyquist_x;
    const int cutoff_x = static_cast<int>(nx * crescent::constants::ONE_THIRD_CUTOFF);
    const int cutoff_y = static_cast<int>(ny * crescent::constants::ONE_THIRD_CUTOFF);
    const Real normalization_factor = Real(nx * ny);

    const crescent::spectral::filter_kernels::CombinedSpectralFilterParameters
        filter_params{
            nx,
            ny,
            dk_x,
            dk_y,
            normalization_factor,
            k_nyquist,
            cutoff_x,
            cutoff_y,
            k_peak,
            beta1,
            beta2,
            rolloff_alpha,
            rolloff_power,
            true
        };

    weights.update(
        *fft,
        1,
        {Real(nx), Real(ny), dk_x, dk_y, k_nyquist, Real(cutoff_x),
         Real(cutoff_y), k_peak, beta1, beta2, rolloff_alpha,
         Real(rolloff_power)},
        [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            return crescent::spectral::filter_kernels::combinedFilterWeight(
                i, j, filter_params);
        });
    const auto weight = weights.localWeights();

    const auto spectral_filter =
        [=] AMREX_GPU_DEVICE(int i, int j, int k, GpuComplex<Real>& fft_c) noexcept
        {
            fft_c = fft_c * weight(i, j, k) / normalization_factor;
        };

    crescent::spectral::applyInPlaceSpectralFilter(eta, fft, spectral_filter);
    crescent::spectral::applyInPlaceSpectralFilter(phi, fft, spectral_filter);
}

} // namespace crescent::nested::filter

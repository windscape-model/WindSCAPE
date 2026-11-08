#include "CRESCENT_Spectral.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_SpectralFilterExecution.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_SpectralFilterKernels.H"
#include "CRESCENT_SpectralWorkspace.H"

#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuComplex.H>
#include <cmath>

namespace crescent::spectral {

namespace FFT = amrex::FFT;
namespace Gpu = amrex::Gpu;
namespace ParallelDescriptor = amrex::ParallelDescriptor;

using amrex::Box;
using amrex::cMultiFab;
using amrex::GeometryData;
using amrex::GpuComplex;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;

// Transform sizes come from the box of the transformed MultiFab: a child
// level covers only part of the problem domain of its geometry.

void dealiasing(MultiFab& mf,
                std::unique_ptr<FFT::R2C<Real>>& fft,
                const GeometryData& geom,
                MultiFab& scratch)
{
    const auto spectral_grid = crescent::shared::makeSpectralGridMetrics(mf, geom);
    int nx = spectral_grid.nx;
    int ny = spectral_grid.ny;

    // The ghost cells of the scratch must be zero: the copy below writes them
    // into mf's ghost region. The transform writes every valid cell.
    scratch.setBndry(0.0);

    // Two-thirds rule on the spectral coefficient fft_x of mode (kx, ky).
    fft->forwardThenBackward(mf, scratch,
        [=] AMREX_GPU_DEVICE (int i, int j, int, GpuComplex<Real>& fft_x) noexcept
        {
            const int kx =
                crescent::spectral::filter_kernels::spectralKxIndex(i, nx);
            const int ky =
                crescent::spectral::filter_kernels::spectralKyIndex(j, ny);

            // 2/3 rule cutoffs in both directions
            int cutoff_x = static_cast<int>(
                nx * crescent::constants::ONE_THIRD_CUTOFF);
            int cutoff_y = static_cast<int>(
                ny * crescent::constants::ONE_THIRD_CUTOFF);

            // AMReX forwardThenBackward gives N² factor, must divide by N² for no net scaling
            fft_x = fft_x * crescent::constants::fftNormalization(nx, ny);

            // Zero out high-frequency modes in BOTH x and y directions (2D filter)
            if (crescent::spectral::filter_kernels::spectralModeExceedsCutoff(
                    kx, ky, cutoff_x, cutoff_y)) {
                fft_x = GpuComplex<Real>(0.0, 0.0);
            }
        }
    );

    MultiFab::Copy(mf, scratch, 0, 0, mf.nComp(), mf.nGrow());
}

void forwardFFT(std::unique_ptr<FFT::R2C<Real>>& fft,
                const GeometryData& geomdata,
                const MultiFab& mf_in,
                cMultiFab& mf_out)
{
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(mf_in, geomdata);
    int nx = spectral_grid.nx;
    int ny = spectral_grid.ny;
    const Real nxy = nx * ny;
    fft->forward(mf_in, mf_out);
    mf_out.mult(1.0/nxy, 0, 1);
}

void backwardFFT(std::unique_ptr<FFT::R2C<Real>>& fft,
                 cMultiFab& mf_in,
                 MultiFab& mf_out)
{
    fft->backward(mf_in, mf_out);
}

void computeWavenumberPowers(MultiFab& wavenumber_pow,
                            const Real dk_x,
                            const Real dk_y)
{
    int ny = wavenumber_pow.boxArray()[0].length(1);

    for (MFIter mfi(wavenumber_pow, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const auto &box = mfi.tilebox();
        auto wavenumber_arr = wavenumber_pow.array(mfi);

        ParallelFor(box, wavenumber_pow.nComp(), [=] AMREX_GPU_DEVICE(int i, int j, int z, int k)
        {
            const Real an =
                crescent::spectral::filter_kernels::spectralMagnitude(
                    i,
                    j,
                    ny,
                    dk_x,
                    dk_y);
            wavenumber_arr(i, j, z, k) = std::pow(an, (k+1));
        });
    }
}

void applyLowPassFilter(std::unique_ptr<FFT::R2C<Real>>& fft,
                        const GeometryData& geomdata,
                        MultiFab& eta,
                        MultiFab& phi,
                        const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
                        crescent::shared::SpectralWeightTable& weights)
{
    const Real dk_x = spectral_settings.wavenumbers.dk_x;
    const Real dk_y = spectral_settings.wavenumbers.dk_y;
    const Real beta1 = spectral_settings.filter.beta1_filter;
    const Real beta2 = spectral_settings.filter.beta2_filter;
    const Real peak_wavenumber = spectral_settings.filter.peak_wavenumber_physical;

    // Exponential low-pass filter about the peak wavenumber. The peak is the
    // initializer's carrier or spectral peak, or the explicit
    // crescent.peak_wavenumber_override; both are positive.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        peak_wavenumber > 0.0,
        "applyLowPassFilter requires a positive peak wavenumber");
    const Real kp = peak_wavenumber;

    const auto spectral_grid = crescent::shared::makeSpectralGridMetrics(eta, geomdata);
    int nx = spectral_grid.nx;
    int ny = spectral_grid.ny;
    // require_positive_peak = true, as on the child path: a vanishing peak
    // would make the low-pass weight exp(-inf) = 0 for every nonzero mode.
    const crescent::spectral::filter_kernels::LowPassFilterParameters
        filter_params{
            ny,
            dk_x,
            dk_y,
            Real(nx * ny),
            kp,
            beta1,
            beta2,
            true,
            true
        };

    weights.update(
        *fft,
        1,
        {Real(nx), Real(ny), dk_x, dk_y, kp, beta1, beta2},
        [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            return crescent::spectral::filter_kernels::lowPassModeWeight(
                i, j, filter_params);
        });
    const auto weight = weights.localWeights();
    const Real normalization_factor = filter_params.normalization_factor;

    // The zero mode is removed. Below the peak the weight is 1, so
    // fft_c * 1 / N equals fft_c / N; the positive peak keeps k = 0 below it.
    const auto spectral_filter =
        [=] AMREX_GPU_DEVICE(int i, int j, int k, GpuComplex<Real>& fft_c) noexcept
        {
            if (i == 0 && j == 0) {
                fft_c = GpuComplex<Real>(0.0, 0.0);
            } else {
                fft_c = fft_c * weight(i, j, k) / normalization_factor;
            }
        };

    crescent::spectral::applySurfacePairSpectralFilter(
        eta,
        phi,
        fft,
        spectral_filter);
}

} // namespace crescent::spectral

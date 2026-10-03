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

// NOTE: Removed getFFTGridDimensions() helper function
// For nested AMR, we must use the actual MultiFab box dimensions,
// NOT the FFT grid or geometry domain dimensions!

void dealiasing(MultiFab& mf,
                std::unique_ptr<FFT::R2C<Real>>& fft,
                const GeometryData& geom,
                MultiFab& scratch)
{
    // CRITICAL: For nested AMR, use the ACTUAL MultiFab grid size, not geometry domain!
    // The geometry domain can be larger than the actual refined region.
    const auto spectral_grid = crescent::shared::makeSpectralGridMetrics(mf, geom);
    int nx = spectral_grid.nx;
    int ny = spectral_grid.ny;

    // Match the fresh-allocation semantics of the prior implementation: the
    // scratch ghost cells must be zero so MultiFab::Copy below does not
    // propagate stale data into mf's ghost region.
    scratch.setVal(0.0);

    // Apply 2/3 rule dealiasing in spectral space
    // NOTE: Despite naming convention, fft_x is a 2D spectral coefficient depending on BOTH kx and ky
    fft->forwardThenBackward(mf, scratch,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, GpuComplex<Real>& fft_x) noexcept
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
    // CRITICAL: Use actual MultiFab grid size for normalization
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(mf_in, geomdata);
    int nx = spectral_grid.nx;
    int ny = spectral_grid.ny;
    const Real nxy = nx * ny;
    fft->forward(mf_in, mf_out);
    mf_out.mult(1.0/nxy, 0, 1);
}

void backwardFFT(std::unique_ptr<FFT::R2C<Real>>& fft,
                 const GeometryData& geomdata,
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
        const auto &box = mfi.validbox();
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
                        cMultiFab& eta_spectral,
                        const Real dk_x,
                        const Real dk_y,
                        const Real beta1,
                        const Real beta2,
                        const Real peak_wavenumber,
                        cMultiFab& peak_wavenumber_gather_scratch,
                        crescent::shared::PeakWavenumberReductionScratch& peak_wavenumber_reduction_scratch)
{
    // Apply exponential low-pass filter with automatic or user-specified peak wavenumber
    Real kp = 0.0;

    // Check if peak wavenumber is provided by user
    if (peak_wavenumber >= 0.0) {
        kp = peak_wavenumber;
    } else {
        forwardFFT(fft, geomdata, eta, eta_spectral);

        // Compute peak wavenumber from spectral energy distribution using the
        // pre-allocated rank-0 gather scratch (shares RootLayouts::cba/cdmap).
        cMultiFab& eta_c_root = peak_wavenumber_gather_scratch;
        ParallelCopy(eta_c_root, eta_spectral, 0, 0, 1);

        Real sum_weights = 0.0;
        auto& d_kp = peak_wavenumber_reduction_scratch.device;
        auto& h_kp = peak_wavenumber_reduction_scratch.host;
        h_kp[0] = 0.0;
        h_kp[1] = 0.0;
        Gpu::copy(Gpu::hostToDevice, h_kp.begin(), h_kp.end(), d_kp.begin());
        Real* d_kp_ptr = d_kp.data();

        if (ParallelDescriptor::IOProcessor()) {
            // CRITICAL: Use actual MultiFab grid size for peak wavenumber calculation
            const auto spectral_grid =
                crescent::shared::makeSpectralGridMetrics(eta, geomdata);
            int nx = spectral_grid.nx;
            int ny = spectral_grid.ny;

            for (MFIter mfi(eta_c_root, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.tilebox();
                bx.growHi(0, -1);
                bx.setRange(1, 0, ny / 2);
                auto eta_arr = eta_c_root.const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                {
                    const Real k_mag =
                        crescent::spectral::filter_kernels::spectralMagnitude(
                            i, j, ny, dk_x, dk_y);

                    Real mag_sq = std::pow(eta_arr(i,j,k).real(), 2) + std::pow(eta_arr(i,j,k).imag(), 2);
                    Real spectral_weight = mag_sq * 4.0;

                    Real local_kp_contrib = std::pow(spectral_weight, 4) * k_mag;
                    Real local_sum = std::pow(spectral_weight, 4);
                    Gpu::Atomic::Add(&d_kp_ptr[0], local_kp_contrib);
                    Gpu::Atomic::Add(&d_kp_ptr[1], local_sum);
                });
            }

            amrex::Gpu::streamSynchronize();

            Gpu::copy(Gpu::deviceToHost, d_kp.begin(), d_kp.end(), h_kp.begin());
            kp = h_kp[0];
            sum_weights = h_kp[1];

            if (sum_weights > 0.0) {
                kp = kp / sum_weights;
            } else {
                kp = 0.0;
            }
        }

        ParallelDescriptor::Bcast(&kp, sizeof(Real), ParallelDescriptor::IOProcessorNumber());
    }

    // Apply filter to eta field
    // CRITICAL: Use actual MultiFab grid size for normalization
    const auto spectral_grid = crescent::shared::makeSpectralGridMetrics(eta, geomdata);
    int nx = spectral_grid.nx;
    int ny = spectral_grid.ny;
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
            false
        };

    const auto spectral_filter =
        [=] AMREX_GPU_DEVICE(int i, int j, int k, GpuComplex<Real>& fft_c) noexcept
        {
            crescent::spectral::filter_kernels::applyLowPassZeroMean(
                i,
                j,
                fft_c,
                filter_params);
        };

    crescent::spectral::applySurfacePairSpectralFilter(
        eta,
        phi,
        fft,
        spectral_filter);
}

} // namespace crescent::spectral

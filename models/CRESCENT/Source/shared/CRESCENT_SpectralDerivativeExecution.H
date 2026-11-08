#pragma once

#include <AMReX_FFT.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>

#include <memory>

#include "CRESCENT_SpectralFilterKernels.H"
#include "CRESCENT_SpectralGridMetrics.H"

namespace crescent::spectral {

template <typename SpectralWeightFunction>
void buildDirectionalDerivativeSpectra(
    amrex::cMultiFab& field_spectral,
    amrex::cMultiFab& field_x_spectral,
    amrex::cMultiFab& field_y_spectral,
    const amrex::MultiFab& field,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    amrex::Real dk_x,
    amrex::Real dk_y,
    SpectralWeightFunction spectral_weight)
{
    fft->forward(field, field_spectral);

    const amrex::Real normalization_factor =
        amrex::Real(spectral_grid.nx * spectral_grid.ny);

    for (amrex::MFIter mfi(field_spectral, amrex::TilingIfNotGPU());
         mfi.isValid();
         ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        auto field_spectral_arr = field_spectral.const_array(mfi);
        auto field_x_spectral_arr = field_x_spectral.array(mfi);
        auto field_y_spectral_arr = field_y_spectral.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const int ky_index =
                crescent::spectral::filter_kernels::spectralKyIndex(
                    j,
                    spectral_grid.ny);
            const amrex::Real sigma = spectral_weight(i, ky_index);

            if (sigma == 0.0) {
                field_x_spectral_arr(i, j, k) = amrex::GpuComplex<amrex::Real>(0.0, 0.0);
                field_y_spectral_arr(i, j, k) = amrex::GpuComplex<amrex::Real>(0.0, 0.0);
                return;
            }

            const amrex::Real weighted_norm = sigma / normalization_factor;

            // Preserve the existing R2C half-spectrum convention in x:
            // the nonnegative i index is the active kx mode index here.
            const amrex::GpuComplex<amrex::Real> phase_x(0.0, i * dk_x);
            const amrex::GpuComplex<amrex::Real> phase_y(0.0, ky_index * dk_y);

            field_x_spectral_arr(i, j, k) =
                field_spectral_arr(i, j, k) * phase_x * weighted_norm;
            field_y_spectral_arr(i, j, k) =
                field_spectral_arr(i, j, k) * phase_y * weighted_norm;
        });
    }
}

template <typename SpectralWeightFunction>
void executeDirectionalDerivativePair(
    amrex::MultiFab& field_x_out,
    amrex::MultiFab& field_y_out,
    const amrex::MultiFab& field,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    amrex::Real dk_x,
    amrex::Real dk_y,
    SpectralWeightFunction spectral_weight,
    amrex::cMultiFab& input_spectral_scratch,
    amrex::cMultiFab& x_spectral_scratch,
    amrex::cMultiFab& y_spectral_scratch)
{
    buildDirectionalDerivativeSpectra(
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch,
        field,
        fft,
        spectral_grid,
        dk_x,
        dk_y,
        spectral_weight);

    fft->backward(x_spectral_scratch, field_x_out);
    fft->backward(y_spectral_scratch, field_y_out);
}

template <typename SpectralWeightFunction>
void executeDirectionalDerivativePairPreservingGrowCells(
    amrex::MultiFab& field_x_out,
    amrex::MultiFab& field_y_out,
    const amrex::MultiFab& field,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const crescent::shared::SpectralGridMetrics& spectral_grid,
    amrex::Real dk_x,
    amrex::Real dk_y,
    SpectralWeightFunction spectral_weight,
    amrex::cMultiFab& input_spectral_scratch,
    amrex::cMultiFab& x_spectral_scratch,
    amrex::cMultiFab& y_spectral_scratch,
    amrex::MultiFab& x_valid_scratch,
    amrex::MultiFab& y_valid_scratch)
{
    executeDirectionalDerivativePair(
        x_valid_scratch,
        y_valid_scratch,
        field,
        fft,
        spectral_grid,
        dk_x,
        dk_y,
        spectral_weight,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);

    amrex::MultiFab::Copy(field_x_out, x_valid_scratch, 0, 0, 1, 0);
    amrex::MultiFab::Copy(field_y_out, y_valid_scratch, 0, 0, 1, 0);
}

} // namespace crescent::spectral

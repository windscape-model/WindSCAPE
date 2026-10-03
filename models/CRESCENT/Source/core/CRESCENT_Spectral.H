#pragma once

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#include "CRESCENT_FreeSurfaceViews.H"
#include "CRESCENT_WaveState.H"

namespace crescent::shared {
struct PeakWavenumberReductionScratch;
}

namespace crescent::spectral {

/**
 * @brief Apply 2/3 dealiasing rule to remove high-frequency modes
 */
void dealiasing(amrex::MultiFab& mf,
                std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                const amrex::GeometryData& geom,
                amrex::MultiFab& scratch);

/**
 * @brief Forward FFT transform
 */
void forwardFFT(std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                const amrex::GeometryData& geomdata,
                const amrex::MultiFab& mf_in,
                amrex::cMultiFab& mf_out);

/**
 * @brief Backward FFT transform
 */
void backwardFFT(std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                 const amrex::GeometryData& geomdata,
                 amrex::cMultiFab& mf_in,
                 amrex::MultiFab& mf_out);

/**
 * @brief Compute wavenumber powers |k|^n for perturbation expansion
 */
void computeWavenumberPowers(amrex::MultiFab& wavenumber_pow,
                            const amrex::Real dk_x,
                            const amrex::Real dk_y);

/**
 * @brief Apply low-pass filter with exponential decay
 *
 * Applies exponential filter: filter(k) = exp(-|k/(β₁·k_peak)|^β₂) for k > k_peak
 * Automatically detects peak wavenumber from spectral energy if not provided.
 *
 * @param fft FFT object
 * @param geomdata Geometry data
 * @param eta Surface elevation field (filtered in-place)
 * @param phi Velocity potential field (filtered in-place)
 * @param eta_spectral Temporary workspace for eta FFT when peak detection is needed
 * @param dk_x Wavenumber spacing in x
 * @param dk_y Wavenumber spacing in y
 * @param beta1 Filter width parameter (typically 8.0)
 * @param beta2 Filter order parameter (typically 12.0)
 * @param peak_wavenumber Optional peak wavenumber override (-1.0 = auto-detect)
 */
void applyLowPassFilter(std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                        const amrex::GeometryData& geomdata,
                        amrex::MultiFab& eta,
                        amrex::MultiFab& phi,
                        amrex::cMultiFab& eta_spectral,
                        const amrex::Real dk_x,
                        const amrex::Real dk_y,
                        const amrex::Real beta1,
                        const amrex::Real beta2,
                        const amrex::Real peak_wavenumber,
                        amrex::cMultiFab& peak_wavenumber_gather_scratch,
                        crescent::shared::PeakWavenumberReductionScratch& peak_wavenumber_reduction_scratch);

inline void applyLowPassFilter(
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geomdata,
    amrex::MultiFab& eta,
    amrex::MultiFab& phi,
    amrex::cMultiFab& eta_spectral,
    const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
    amrex::cMultiFab& peak_wavenumber_gather_scratch,
    crescent::shared::PeakWavenumberReductionScratch& peak_wavenumber_reduction_scratch)
{
    applyLowPassFilter(
        fft,
        geomdata,
        eta,
        phi,
        eta_spectral,
        spectral_settings.wavenumbers.dk_x,
        spectral_settings.wavenumbers.dk_y,
        spectral_settings.filter.beta1_filter,
        spectral_settings.filter.beta2_filter,
        spectral_settings.filter.peak_wavenumber_physical,
        peak_wavenumber_gather_scratch,
        peak_wavenumber_reduction_scratch);
}

inline void applyLowPassFilter(
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geomdata,
    crescent::shared::FreeSurfaceStateView state,
    amrex::cMultiFab& eta_spectral,
    const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
    amrex::cMultiFab& peak_wavenumber_gather_scratch,
    crescent::shared::PeakWavenumberReductionScratch& peak_wavenumber_reduction_scratch)
{
    applyLowPassFilter(
        fft,
        geomdata,
        state.eta,
        state.phi,
        eta_spectral,
        spectral_settings,
        peak_wavenumber_gather_scratch,
        peak_wavenumber_reduction_scratch);
}

} // namespace crescent::spectral

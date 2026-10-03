#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_FFT.H>
#include "CRESCENT_Constants.H"
#include "CRESCENT_WaveState.H"

namespace crescent::nested::filter {

/**
 * @brief Apply Hou-Li exponential filter to L1 fields (eta, phi) in spectral space
 *
 * Purpose: Remove high-frequency modes from L1 state variables after RK4 update
 *          to prevent energy accumulation and maintain stability.
 *
 * This function:
 * 1. Forward FFT: eta, phi → spectral space
 * 2. Apply combined filter:
 *    - Hou-Li exponential filter (grid-scale stabilization)
 *    - 2/3 dealiasing rule (zero modes beyond 2/3 k_Nyquist)
 *    - Optional: low-pass filter (zero modes beyond k_peak)
 * 3. Backward FFT: spectral → physical space
 *
 * @param eta        Surface elevation field (modified in-place)
 * @param phi        Velocity potential field (modified in-place)
 * @param fft        FFT object for L1 grid
 * @param geomdata   L1 geometry data
 * @param dk_x       Wavenumber spacing in x (rad/m)
 * @param dk_y       Wavenumber spacing in y (rad/m)
 * @param k_peak     Physical peak wavenumber (rad/m), -1 to disable low-pass
 * @param beta1      Low-pass filter parameter (default: 9.0)
 * @param beta2      Low-pass filter exponent (default: 12.0)
 * @param log_diagnostics Emit one caller-selected diagnostic block
 */
void applyHouLiFilter_L1(
    amrex::MultiFab& eta,
    amrex::MultiFab& phi,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geomdata,
    amrex::Real dk_x,
    amrex::Real dk_y,
    amrex::Real k_peak = -1.0,
    amrex::Real beta1 = 9.0,
    amrex::Real beta2 = 12.0,
    bool log_diagnostics = false
);

inline void applyHouLiFilter_L1(
    amrex::MultiFab& eta,
    amrex::MultiFab& phi,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geomdata,
    const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
    bool log_diagnostics = false)
{
    applyHouLiFilter_L1(
        eta,
        phi,
        fft,
        geomdata,
        spectral_settings.wavenumbers.dk_x,
        spectral_settings.wavenumbers.dk_y,
        spectral_settings.filter.peak_wavenumber_physical,
        spectral_settings.filter.beta1_filter,
        spectral_settings.filter.beta2_filter,
        log_diagnostics);
}

} // namespace crescent::nested::filter

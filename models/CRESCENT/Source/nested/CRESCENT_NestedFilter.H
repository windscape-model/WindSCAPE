#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_FFT.H>
#include "CRESCENT_SpectralWeightTable.H"
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
 * @param spectral_settings Child wavenumber spacing, peak, low-pass and roll-off parameters
 * @param weights    Cached per-mode filter weights of the level
 *
 * The k = 0 mode is kept: on the child window neither eta nor phi has zero
 * mean, and the state blend mixes the filtered copy against interpolated
 * parent data with nonzero mean.
 */
void applyHouLiFilter_L1(
    amrex::MultiFab& eta,
    amrex::MultiFab& phi,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geomdata,
    const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
    crescent::shared::SpectralWeightTable& weights);

} // namespace crescent::nested::filter

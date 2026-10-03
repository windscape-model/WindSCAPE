#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>

#include "CRESCENT_FreeSurfaceViews.H"

namespace crescent::init {

/**
 * @brief Initialize Stokes wave field
 *
 * Creates a 2nd-order Stokes wave on the root processor and distributes to all processors.
 *
 * @param eta Output: surface elevation field (distributed)
 * @param phi Output: velocity potential field (distributed)
 * @param ba_root Root box array (rank 0 only)
 * @param dmap_root Root distribution mapping (rank 0 only)
 * @param geom Geometry data
 * @param dk_x Wavenumber spacing in x
 * @param dk_y Wavenumber spacing in y
 * @param wave_steepness Wave steepness parameter
 * @param inverse_gravity Inverse gravity (1/g)
 * @param wave_mode_x Wave mode index - sets total wavenumber magnitude
 * @param wave_angle Wave propagation angle in radians (0 = +x, π/2 = +y)
 */
void initializeStokesWave(amrex::MultiFab& eta,
                          amrex::MultiFab& phi,
                          const amrex::BoxArray& ba_root,
                          const amrex::DistributionMapping& dmap_root,
                          const amrex::GeometryData& geom,
                          const amrex::Real dk_x,
                          const amrex::Real dk_y,
                          const amrex::Real wave_steepness,
                          const amrex::Real inverse_gravity,
                          const int wave_mode_x,
                          const amrex::Real wave_angle);

inline void initializeStokesWave(
    crescent::shared::FreeSurfaceStateView state,
    const amrex::BoxArray& ba_root,
    const amrex::DistributionMapping& dmap_root,
    const amrex::GeometryData& geom,
    const amrex::Real dk_x,
    const amrex::Real dk_y,
    const amrex::Real wave_steepness,
    const amrex::Real inverse_gravity,
    const int wave_mode_x,
    const amrex::Real wave_angle)
{
    initializeStokesWave(
        state.eta,
        state.phi,
        ba_root,
        dmap_root,
        geom,
        dk_x,
        dk_y,
        wave_steepness,
        inverse_gravity,
        wave_mode_x,
        wave_angle);
}

/**
 * @brief Generate JONSWAP spectrum waves
 *
 * Generates a wave field following the JONSWAP spectrum with directional spreading.
 * Initializes in spectral space, transforms to physical space, and applies dealiasing.
 *
 * @param eta_saved Output: surface elevation field
 * @param phi_saved Output: velocity potential field
 * @param fft FFT object for transforms
 * @param geom Geometry data
 * @param cba_root Complex box array for root processor
 * @param cdmap_root Complex distribution mapping for root processor
 * @param dk_x Wavenumber spacing in x
 * @param dk_y Wavenumber spacing in y
 * @param wind_speed_10m Wind speed at 10m height
 * @param fetch Fetch length
 * @param inverse_gravity Inverse gravity (1/g)
 * @param initial_angle Initial wave propagation angle
 * @param angle_spread Angular spreading parameter
 * @param directional_exponent Exponent for directional spreading
 * @param alpha Output: JONSWAP alpha parameter (computed)
 * @param peak_frequency Output: peak frequency (computed)
 * @param peak_wavelength Output: peak wavelength (computed)
 */
void generateJonswapWaves(amrex::MultiFab& eta_saved,
                          amrex::MultiFab& phi_saved,
                          std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                          const amrex::GeometryData& geom,
                          const amrex::BoxArray& cba_root,
                          const amrex::DistributionMapping& cdmap_root,
                          const amrex::Real dk_x,
                          const amrex::Real dk_y,
                          const amrex::Real wind_speed_10m,
                          const amrex::Real fetch,
                          const amrex::Real inverse_gravity,
                          const amrex::Real initial_angle,
                          const amrex::Real angle_spread,
                          const int directional_exponent,
                          amrex::Real& alpha,
                          amrex::Real& peak_frequency,
                          amrex::Real& peak_wavelength,
                          int verbosity);

inline void generateJonswapWaves(
    crescent::shared::FreeSurfaceStateView state,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geom,
    const amrex::BoxArray& cba_root,
    const amrex::DistributionMapping& cdmap_root,
    const amrex::Real dk_x,
    const amrex::Real dk_y,
    const amrex::Real wind_speed_10m,
    const amrex::Real fetch,
    const amrex::Real inverse_gravity,
    const amrex::Real initial_angle,
    const amrex::Real angle_spread,
    const int directional_exponent,
    amrex::Real& alpha,
    amrex::Real& peak_frequency,
    amrex::Real& peak_wavelength,
    int verbosity)
{
    generateJonswapWaves(
        state.eta,
        state.phi,
        fft,
        geom,
        cba_root,
        cdmap_root,
        dk_x,
        dk_y,
        wind_speed_10m,
        fetch,
        inverse_gravity,
        initial_angle,
        angle_spread,
        directional_exponent,
        alpha,
        peak_frequency,
        peak_wavelength,
        verbosity);
}

} // namespace crescent::init

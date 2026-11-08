#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>


namespace crescent::init {

/**
 * @brief Initialize Stokes wave field
 *
 * Creates a linear (Airy) progressive wave on the root processor and
 * distributes it to all processors.
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
 * @param wave_mode_x Absolute x-mode away from the y axis
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

/**
 * @brief Initialize a crescent (3D horseshoe) wave field
 *
 * Writes a 5th-order Stokes carrier along +x plus an oblique free-wave
 * perturbation (the class-II instability seed) on the root processor and
 * distributes to all processors. Ported verbatim (signs included) from the
 * reference Fortran HOS crescent_init, generalized to arbitrary gravity.
 *
 * @param eta Output: surface elevation field (distributed)
 * @param phi Output: velocity potential field (distributed)
 * @param ba_root Root box array (rank 0 only)
 * @param dmap_root Root distribution mapping (rank 0 only)
 * @param geom Geometry data
 * @param dk_x Wavenumber spacing in x
 * @param dk_y Wavenumber spacing in y
 * @param wave_steepness Carrier steepness ak
 * @param inverse_gravity Inverse gravity (1/g)
 * @param wave_mode_x Carrier x-mode index (carrier k = wave_mode_x * dk_x)
 * @param crescent_eps Perturbation amplitude relative to carrier amplitude
 * @param crescent_kx_ratio Perturbation kx / carrier k
 * @param crescent_ky_ratio Perturbation ky / carrier k
 * @param crescent_stokes_order Stokes order of the analytic carrier (1..5)
 */
void initializeCrescentWave(amrex::MultiFab& eta,
                            amrex::MultiFab& phi,
                            const amrex::BoxArray& ba_root,
                            const amrex::DistributionMapping& dmap_root,
                            const amrex::GeometryData& geom,
                            const amrex::Real dk_x,
                            const amrex::Real dk_y,
                            const amrex::Real wave_steepness,
                            const amrex::Real inverse_gravity,
                            const int wave_mode_x,
                            const amrex::Real crescent_eps,
                            const amrex::Real crescent_kx_ratio,
                            const amrex::Real crescent_ky_ratio,
                            const int crescent_stokes_order);

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
 * @param use_directional_spectrum Whether to use the directional spectrum
 * @param jonswap_amplitude_scale Multiplier on spectrum-derived Fourier amplitudes
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
                          const int use_directional_spectrum,
                          const amrex::Real jonswap_amplitude_scale,
                          amrex::Real& alpha,
                          amrex::Real& peak_frequency,
                          amrex::Real& peak_wavelength);

} // namespace crescent::init

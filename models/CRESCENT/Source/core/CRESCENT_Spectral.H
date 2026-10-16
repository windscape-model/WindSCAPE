#pragma once

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#include "CRESCENT_SpectralWeightTable.H"
#include "CRESCENT_WaveState.H"

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
                 amrex::cMultiFab& mf_in,
                 amrex::MultiFab& mf_out);

/**
 * @brief Compute wavenumber powers |k|^n for perturbation expansion
 */
void computeWavenumberPowers(amrex::MultiFab& wavenumber_pow,
                            const amrex::Real dk_x,
                            const amrex::Real dk_y);

/**
 * @brief Exponential low-pass filter of a surface pair, in place
 *
 * Weight exp(-|k/(beta1 k_peak)|^beta2) above the peak wavenumber k_peak,
 * with the zero mode removed (the parent's filter and dealiasing pass).
 */
void applyLowPassFilter(std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                        const amrex::GeometryData& geomdata,
                        amrex::MultiFab& eta,
                        amrex::MultiFab& phi,
                        const crescent::shared::LevelWaveSpectralSettings& spectral_settings,
                        crescent::shared::SpectralWeightTable& weights);

} // namespace crescent::spectral

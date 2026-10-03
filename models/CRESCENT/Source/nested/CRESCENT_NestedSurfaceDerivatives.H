#ifndef CRESCENT_NESTED_SURFACE_DERIVATIVES_H
#define CRESCENT_NESTED_SURFACE_DERIVATIVES_H

#include <AMReX_Box.H>
#include <AMReX_FFT_R2C.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace crescent::nested::derivatives {

/**
 * @brief Compute periodic L1 free-surface gradients on the full nested grid
 *
 * Algorithm:
 *   1. Operate on the prepared full L1 periodic copies
 *   2. FFT with L1 dimensions and L1 wavenumber spacing
 *   3. Compute x/y surface derivatives
 *   4. Return full L1 derivative fields (Ω_S extraction remains a caller concern)
 */
void computeSurfaceDerivatives(
    amrex::MultiFab& eta_x,
    amrex::MultiFab& eta_y,
    amrex::MultiFab& phi_x,
    amrex::MultiFab& phi_y,
    const amrex::MultiFab& eta,
    const amrex::MultiFab& phi,
    const amrex::Box& box_S,
    const amrex::Geometry& geom,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    int spectral_filter_mode,
    amrex::Real k_peak,
    amrex::Real beta1,
    amrex::Real beta2,
    amrex::cMultiFab& input_spectral_scratch,
    amrex::cMultiFab& x_spectral_scratch,
    amrex::cMultiFab& y_spectral_scratch,
    amrex::MultiFab& x_valid_scratch,
    amrex::MultiFab& y_valid_scratch,
    int verbosity);

} // namespace crescent::nested::derivatives

#endif // CRESCENT_NESTED_SURFACE_DERIVATIVES_H

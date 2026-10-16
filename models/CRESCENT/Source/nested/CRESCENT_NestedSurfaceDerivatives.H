#ifndef CRESCENT_NESTED_SURFACE_DERIVATIVES_H
#define CRESCENT_NESTED_SURFACE_DERIVATIVES_H

#include <AMReX_FFT_R2C.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace crescent::nested::derivatives {

/**
 * @brief Spectral x/y derivatives of the periodic child copies.
 *
 * One forward transform per field, the square two-thirds cut applied to the
 * derivative spectra, two inverse transforms. The wavenumber spacing is
 * that of the child box (n dx), not of the parent domain.
 */
void computeSurfaceDerivatives(
    amrex::MultiFab& eta_x,
    amrex::MultiFab& eta_y,
    amrex::MultiFab& phi_x,
    amrex::MultiFab& phi_y,
    const amrex::MultiFab& eta,
    const amrex::MultiFab& phi,
    const amrex::Geometry& geom,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    amrex::cMultiFab& input_spectral_scratch,
    amrex::cMultiFab& x_spectral_scratch,
    amrex::cMultiFab& y_spectral_scratch,
    amrex::MultiFab& x_valid_scratch,
    amrex::MultiFab& y_valid_scratch);

} // namespace crescent::nested::derivatives

#endif // CRESCENT_NESTED_SURFACE_DERIVATIVES_H

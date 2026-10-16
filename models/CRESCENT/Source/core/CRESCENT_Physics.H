#pragma once

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_SpectralWeightTable.H"

namespace crescent::physics {

/**
 * @brief Surface slopes and potential gradients from one forward FFT per field
 *
 * Computes eta_x, eta_y, phi_x, phi_y on a periodic level.
 */
void computeSurfaceDerivatives(amrex::MultiFab& eta_x, amrex::MultiFab& eta_y,
                                amrex::MultiFab& phi_x, amrex::MultiFab& phi_y,
                                const amrex::MultiFab& eta,
                                const amrex::MultiFab& phi,
                                std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                                const amrex::GeometryData& geom,
                                const amrex::Real dk_x, const amrex::Real dk_y,
                                amrex::cMultiFab& input_spectral_scratch,
                                amrex::cMultiFab& x_spectral_scratch,
                                amrex::cMultiFab& y_spectral_scratch);

/**
 * @brief Scratch arrays of the vertical-velocity series (non-owning)
 */
struct VerticalVelocityWorkspace
{
    amrex::MultiFab& eta_pow;                 // eta^(n+1)/(n+1)!, n = 0..M-2
    amrex::MultiFab& phi_expansion;           // potential coefficients phi_m, m = 0..M-1
    amrex::MultiFab& expansion_deriv;
    amrex::cMultiFab& deriv_spectral;
    amrex::MultiFab& expansion_accum;         // cumulative sums of phi_m
    amrex::MultiFab& spectral_work;
    amrex::cMultiFab& phi_expansion_spectrum;
    amrex::MultiFab& dealiasing_scratch;
    crescent::shared::SpectralWeightTable& kpow_table; // |k|^n, n = 1..M
};

/**
 * @brief Surface vertical velocity phi_z of a periodic pair (eta, phi)
 *
 * Order-M HOS expansion of the Dirichlet-to-Neumann operator (West et al.
 * 1987). wavenumber_pow holds |k|^(n+1) on the spectral layout; dk_x, dk_y
 * are the spectral spacing of the transform box. The final 2/3 projection
 * removes the product modes of the completed sum.
 */
void computeVerticalVelocity(amrex::MultiFab& phi_z,
                              const amrex::MultiFab& eta,
                              const amrex::MultiFab& phi,
                              const amrex::MultiFab& wavenumber_pow,
                              const VerticalVelocityWorkspace& workspace,
                              std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                              const amrex::GeometryData& geom,
                              int perturbation_order,
                              amrex::Real dk_x,
                              amrex::Real dk_y,
                              bool apply_final_projection);

/**
 * @brief Horizontal surface velocities u = phi_x - phi_z eta_x, v = phi_y - phi_z eta_y
 *
 * apply_spectral_dealiasing must be false for blended, non-periodic child
 * fields.
 */
void computeHorizontalVelocities(
    amrex::MultiFab& u,
    amrex::MultiFab& v,
    crescent::shared::DerivativeBundleView derivatives,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geom,
    amrex::MultiFab& dealiasing_scratch,
    bool apply_spectral_dealiasing = true);

} // namespace crescent::physics

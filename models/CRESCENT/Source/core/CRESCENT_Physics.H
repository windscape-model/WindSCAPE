#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_DNOConfig.H"
#include "CRESCENT_DNOWorkspaceViews.H"

namespace crescent::physics {

/**
 * @brief Compute powers of surface elevation for Taylor expansion
 * Computes eta_pow[k] = η^(k+1) / (k+1)! for k = 0 to perturbation_order-2
 *
 * @param eta_pow Output MultiFab with perturbation_order-1 components
 * @param eta Surface elevation field
 * @param fft FFT object for dealiasing
 * @param geom Geometry data
 * @param perturbation_order Order of HOS expansion
 */
void computeEtaPowers(amrex::MultiFab& eta_pow,
                      const amrex::MultiFab& eta,
                      std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                      const amrex::GeometryData& geom,
                      const int perturbation_order,
                      amrex::MultiFab& dealiasing_scratch);

/**
 * @brief Compute surface derivatives using FFT kernel fusion
 * Computes ∂η/∂x, ∂η/∂y, ∂φ/∂x, ∂φ/∂y in a single fused kernel
 *
 * @param eta_x Output: ∂η/∂x
 * @param eta_y Output: ∂η/∂y
 * @param phi_x Output: ∂φ/∂x
 * @param phi_y Output: ∂φ/∂y
 * @param eta Surface elevation
 * @param phi Velocity potential at surface
 * @param fft FFT object
 * @param geom Geometry data
 * @param dk_x Wavenumber spacing in x
 * @param dk_y Wavenumber spacing in y
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
 * @brief Expand velocity potential φ(x,y,z=η) using Taylor expansion
 * Computes φ_k for k = 0 to perturbation_order-1
 *
 * @param phi_expansion Output MultiFab with perturbation_order components
 * @param phi Velocity potential at surface
 * @param wavenumber_pow Wavenumber powers |k|^n
 * @param workspace Non-owning view of the phi-expansion stage workspace
 * @param fft FFT object
 * @param geom Geometry data
 * @param perturbation_order Order of HOS expansion
 */
void expandVelocityPotential(const amrex::MultiFab& phi,
                              const amrex::MultiFab& wavenumber_pow,
                              crescent::dno::PhiExpansionWorkspaceView workspace,
                              std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                              const amrex::GeometryData& geom,
                              const crescent::config::DNOConfig& config);

/**
 * @brief Compute vertical velocity φ_z at surface using HOS expansion
 *
 * @param phi_z Output: vertical velocity at surface
 * @param wavenumber_pow Wavenumber powers
 * @param workspace Non-owning view of the vertical-velocity stage workspace
 * @param fft FFT object
 * @param geom Geometry data
 * @param config DNO configuration (perturbation order, verbosity, etc.)
 */
void computeVerticalVelocity(amrex::MultiFab& phi_z,
                              const amrex::MultiFab& wavenumber_pow,
                              crescent::dno::VerticalVelocityWorkspaceView workspace,
                              std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                              const amrex::GeometryData& geom,
                              const crescent::config::DNOConfig& config,
                              amrex::MultiFab& dealiasing_scratch);

/**
 * @brief Compute horizontal velocities (u, v) at surface
 * u = φ_x - φ_z * η_x
 * v = φ_y - φ_z * η_y
 *
 * @param u Output: horizontal velocity in x
 * @param v Output: horizontal velocity in y
 * @param phi_z Vertical velocity at surface
 * @param eta_x ∂η/∂x
 * @param eta_y ∂η/∂y
 * @param phi_x ∂φ/∂x
 * @param phi_y ∂φ/∂y
 * @param fft FFT object for dealiasing
 * @param geom Geometry data
 */
void computeHorizontalVelocities(amrex::MultiFab& u, amrex::MultiFab& v,
                                  const amrex::MultiFab& phi_z,
                                  const amrex::MultiFab& eta_x, const amrex::MultiFab& eta_y,
                                  const amrex::MultiFab& phi_x, const amrex::MultiFab& phi_y,
                                  std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                                  const amrex::GeometryData& geom,
                                  amrex::MultiFab& dealiasing_scratch);

void computeHorizontalVelocities(
    amrex::MultiFab& u,
    amrex::MultiFab& v,
    crescent::shared::DerivativeBundleView derivatives,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geom,
    amrex::MultiFab& dealiasing_scratch);

} // namespace crescent::physics

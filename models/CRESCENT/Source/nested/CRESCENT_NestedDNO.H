#ifndef CRESCENT_NESTED_DNO_H
#define CRESCENT_NESTED_DNO_H

#include <AMReX_FFT_R2C.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_DNOConfig.H"
#include "CRESCENT_DNOWorkspaceViews.H"

namespace crescent::nested::derivatives {

/**
 * @brief Expand the prepared L1 surface potential through the HOS recurrence
 *
 * Follows the L0 phi-expansion interface, but uses the nested L1 spectral
 * geometry and workspace views.
 */
void expandVelocityPotential(
    const amrex::MultiFab& phi,
    const amrex::MultiFab& wavenumber_pow,
    crescent::dno::PhiExpansionWorkspaceView workspace,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geom,
    const crescent::config::DNOConfig& config);

/**
 * @brief Compute the prepared L1 vertical surface velocity through the HOS DNO
 *
 * Follows the L0 vertical-velocity interface while preserving the nested path's
 * no-intermediate-dealiasing behavior.
 */
void computeVerticalVelocity(
    amrex::MultiFab& phi_z,
    const amrex::MultiFab& wavenumber_pow,
    crescent::dno::VerticalVelocityWorkspaceView workspace,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geom,
    const crescent::config::DNOConfig& config,
    amrex::MultiFab& dealiasing_scratch);

} // namespace crescent::nested::derivatives

#endif // CRESCENT_NESTED_DNO_H

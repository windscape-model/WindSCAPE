#ifndef CRESCENT_NESTED_ENERGY_H
#define CRESCENT_NESTED_ENERGY_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <string>

namespace crescent::nested::energy {

/**
 * @brief Compute total wave energy (potential + kinetic)
 *
 * Energy = ∫∫ (g·η²/2 + ψ·G(η)ψ/2) dx dy, where
 * G(η)ψ = (1 + |grad η|²) w - grad η · grad ψ and w is the surface
 * vertical velocity.
 *
 * @param eta Surface elevation
 * @param psi Velocity potential at surface
 * @param eta_x x derivative of surface elevation
 * @param eta_y y derivative of surface elevation
 * @param psi_x x derivative of surface potential
 * @param psi_y y derivative of surface potential
 * @param w Vertical velocity at surface (not the DNO itself)
 * @param geom Geometry for cell volumes
 * @param gravity Gravitational acceleration (default 9.81)
 * @return Total energy
 */
amrex::Real computeTotalEnergy(
    const amrex::MultiFab& eta,
    const amrex::MultiFab& psi,
    const amrex::MultiFab& eta_x,
    const amrex::MultiFab& eta_y,
    const amrex::MultiFab& psi_x,
    const amrex::MultiFab& psi_y,
    const amrex::MultiFab& w,
    const amrex::Geometry& geom,
    amrex::Real gravity = 9.81
);

/**
 * @brief Compute total wave energy restricted to an index-space region
 *
 * Same integrand as computeTotalEnergy, but only cells inside `region`
 * (in the level's index space) contribute. Used for zone-resolved energy
 * diagnostics (e.g. E over Omega_A only).
 */
amrex::Real computeTotalEnergyInBox(
    const amrex::MultiFab& eta,
    const amrex::MultiFab& psi,
    const amrex::MultiFab& eta_x,
    const amrex::MultiFab& eta_y,
    const amrex::MultiFab& psi_x,
    const amrex::MultiFab& psi_y,
    const amrex::MultiFab& w,
    const amrex::Geometry& geom,
    const amrex::Box& region,
    amrex::Real gravity = 9.81
);

/**
 * @brief Write energy tracking data to CSV file
 *
 * Appends a row with: timestep, time, energy_L0, energy_L1, energy_total,
 * energy_change_percent, energy_L1_omegaA, energy_L0_under_omegaA,
 * energy_composite. Creates file with header if it doesn't exist.
 *
 * energy_composite = energy_L0 - energy_L0_under_omegaA + energy_L1_omegaA
 * is the non-double-counting hierarchy energy: L0 everywhere except the
 * Omega_A footprint, where the refined L1 solution replaces it.
 */
void writeEnergyToCSV(
    const std::string& filename,
    int istep,
    amrex::Real time,
    amrex::Real energy_L0,
    amrex::Real energy_L1,
    amrex::Real energy_initial_total,
    amrex::Real energy_L1_omega_a,
    amrex::Real energy_L0_under_omega_a,
    amrex::Real energy_composite
);

} // namespace crescent::nested::energy

#endif // CRESCENT_NESTED_ENERGY_H

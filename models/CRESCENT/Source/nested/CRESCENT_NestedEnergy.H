#ifndef CRESCENT_NESTED_ENERGY_H
#define CRESCENT_NESTED_ENERGY_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <string>

namespace crescent::nested::energy {

/**
 * @brief Compute total wave energy (potential + kinetic)
 *
 * Energy = ∫∫ (g·η²/2 + ψ·w/2) dx dy
 *
 * @param eta Surface elevation
 * @param psi Velocity potential at surface
 * @param w Vertical velocity at surface
 * @param geom Geometry for cell volumes
 * @param gravity Gravitational acceleration (default 9.81)
 * @return Total energy
 */
amrex::Real computeTotalEnergy(
    const amrex::MultiFab& eta,
    const amrex::MultiFab& psi,
    const amrex::MultiFab& w,
    const amrex::Geometry& geom,
    amrex::Real gravity = 9.81
);

/**
 * @brief Write energy tracking data to CSV file
 *
 * Appends a row with: timestep, time, energy_L0, energy_L1, energy_total
 * Creates file with header if it doesn't exist
 *
 * @param filename Output CSV filename
 * @param istep Timestep number
 * @param time Physical time
 * @param energy_L0 Energy on L0 grid
 * @param energy_L1 Energy on L1 grid (0 if no L1)
 */
void writeEnergyToCSV(
    const std::string& filename,
    int istep,
    amrex::Real time,
    amrex::Real energy_L0,
    amrex::Real energy_L1,
    amrex::Real energy_initial_total
);

} // namespace crescent::nested::energy

#endif // CRESCENT_NESTED_ENERGY_H

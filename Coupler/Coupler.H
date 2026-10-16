#ifndef COUPLER_H
#define COUPLER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>

// No CRESCENT forward declaration needed - use void* for opaque pointer

/**
 * @file Coupler.H
 * @brief ERF-CRESCENT in-process coupling interface for air-sea interaction
 *
 * This header defines the coupling functions that connect ERF atmospheric
 * solver with CRESCENT wave solver for bidirectional exchange of:
 *   - ERF → CRESCENT: Surface atmospheric pressure
 *   - CRESCENT → ERF: Wave elevation and velocities
 *
 * Architecture: In-process coupling (single executable, sequential execution)
 * All MPI ranks run ERF, then all run CRESCENT, maintaining 100% resource utilization.
 *
 * @see the top-level README.md ("Running a coupled case") for the coupling design, and
 *      models/ERF/UPSTREAM.md for the inventory of ERF-side hooks this interface drives.
 */

namespace coupler {


/**
 * @brief Initialize CRESCENT wave solver with parameters from ERF
 *
 * Creates CRESCENT instance, sets up grid matching ERF horizontal domain
 * but with finer resolution (grid_ratio), and initializes wave physics.
 *
 * @param crescent_solver Opaque pointer to CRESCENT instance (will be allocated)
 * @param erf_geom ERF geometry at level 0
 * @param erf_grids ERF grids at level 0
 * @param use_crescent_coupling Output: whether coupling is enabled
 * @param water_density Output: water density for pressure conversion [kg/m³]
 * @param crescent_dt_ratio Output: time step ratio (dt_hos = dt_erf / ratio)
 * @param crescent_grid_ratio Output: grid ratio (dx_hos = dx_erf / ratio)
 */
void InitializeCRESCENTCoupling(
    void*& crescent_solver,
    const amrex::Geometry& erf_geom,
    const amrex::BoxArray& erf_grids,
    bool& use_crescent_coupling,
    amrex::Real& water_density,
    int& crescent_dt_ratio,
    int& crescent_grid_ratio);

/**
 * @brief Transfer atmospheric surface pressure from ERF to CRESCENT
 *
 * Takes surface pressure field, interpolates from coarse ERF grid to fine CRESCENT grid,
 * divides by water density, and passes to CRESCENT for wave evolution.
 *
 * @param crescent_solver Opaque CRESCENT instance pointer
 * @param erf_geom ERF geometry
 * @param erf_grids ERF grids
 * @param erf_dmap ERF distribution mapping
 * @param crescent_grid_ratio Grid refinement ratio
 * @param erf_p_surface Surface pressure field (2D, already extracted at k=0) [Pa]
 * @param time Current simulation time [s]
 * @param lev AMR level (typically 0)
 */
void TransferPressureToCRESCENT(
    void* crescent_solver,
    const amrex::Geometry& erf_geom,
    const amrex::BoxArray& erf_grids,
    const amrex::DistributionMapping& erf_dmap,
    int crescent_grid_ratio,
    const amrex::MultiFab& erf_p_surface,
    amrex::Real time,
    int crescent_lev = 0);

/**
 * @brief Advance CRESCENT with sub-cycling
 *
 * CRESCENT typically requires smaller timestep than ERF for stability.
 * This function advances CRESCENT multiple times (crescent_dt_ratio substeps)
 * to catch up with ERF time.
 *
 * @param crescent_solver Opaque CRESCENT instance pointer
 * @param dt_erf ERF timestep size [s]
 * @param crescent_dt_ratio Number of CRESCENT substeps per ERF step
 */
void AdvanceCRESCENTWithSubcycling(
    void* crescent_solver,
    amrex::Real dt_erf,
    int crescent_dt_ratio);

/**
 * @brief Transfer wave fields from CRESCENT to ERF
 *
 * Extracts surface elevation (η) and velocities (u, v, w) from CRESCENT,
 * coarsens from fine CRESCENT grid to coarse ERF grid, and returns
 * MultiFabs on ERF grid for boundary condition application.
 *
 * @param crescent_solver Opaque CRESCENT instance pointer
 * @param erf_geom ERF geometry
 * @param erf_grids ERF grids
 * @param erf_dmap ERF distribution mapping
 * @param crescent_grid_ratio Grid refinement ratio
 * @param erf_eta Output: Surface elevation on ERF grid [m]
 * @param erf_u Output: Surface u-velocity on ERF grid [m/s]
 * @param erf_v Output: Surface v-velocity on ERF grid [m/s]
 * @param erf_w Output: Surface w-velocity on ERF grid [m/s]
 * @param lev AMR level (typically 0)
 */
void TransferWavesFromCRESCENT(
    void* crescent_solver,
    const amrex::Geometry& erf_geom,
    const amrex::BoxArray& erf_grids,
    const amrex::DistributionMapping& erf_dmap,
    int crescent_grid_ratio,
    amrex::MultiFab& erf_eta,
    amrex::MultiFab& erf_u,
    amrex::MultiFab& erf_v,
    amrex::MultiFab& erf_w,
    int crescent_lev = 0);

/**
 * @brief Destroy CRESCENT solver instance
 *
 * Deallocates CRESCENT solver created by InitializeCRESCENTCoupling.
 * Must be called to avoid memory leak.
 *
 * @param crescent_solver Opaque CRESCENT instance pointer (will be set to nullptr)
 */
void DestroyCRESCENTSolver(void*& crescent_solver);

/**
 * @brief Apply wave boundary condition to ERF terrain
 *
 * Updates ERF's moving terrain (z_phys_nd) with wave elevation from CRESCENT
 * and sets surface velocity boundary conditions. This makes ERF's
 * computational domain follow the wave surface.
 *
 * @param erf_eta Surface elevation [m]
 * @param erf_u Surface u-velocity [m/s]
 * @param erf_v Surface v-velocity [m/s]
 * @param erf_w Surface w-velocity [m/s]
 * @param z_phys_nd ERF terrain height field (nodal)
 * @param erf_geom ERF geometry
 * @param lev AMR level (typically 0)
 */
void ApplyWaveBoundaryCondition(
    const amrex::MultiFab& erf_eta,
    const amrex::MultiFab& erf_u,
    const amrex::MultiFab& erf_v,
    const amrex::MultiFab& erf_w,
    amrex::MultiFab& z_phys_nd,
    const amrex::Geometry& erf_geom,
    int lev = 0);



} // namespace coupler

#endif // COUPLER_H

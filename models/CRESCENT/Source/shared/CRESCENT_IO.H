#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <string>

// Forward declarations
namespace crescent::nested::domains {
    struct DomainStructure;
}

namespace crescent::io {

struct TecplotDerivativeLevels
{
    const amrex::Vector<amrex::MultiFab>& eta_x;
    const amrex::Vector<amrex::MultiFab>& eta_y;
    const amrex::Vector<amrex::MultiFab>& phi_x;
    const amrex::Vector<amrex::MultiFab>& phi_y;
    const amrex::Vector<amrex::MultiFab>& phi_z;
};

struct NestedTecplotDiagnosticLevels
{
    TecplotDerivativeLevels core_surface_derivatives;
    TecplotDerivativeLevels nested_output_derivatives;
    const amrex::Vector<amrex::MultiFab>& chi_D_mask;
    const amrex::Vector<amrex::MultiFab>& sigma_blend;
    const amrex::Vector<amrex::MultiFab>& eta_x_core_reference_l1;
    const amrex::Vector<amrex::MultiFab>& phi_x_core_reference_l1;
};

/**
 * @brief Write AMReX plotfile with wave field data
 *
 * @param plotfilename Base name for the plotfile
 * @param step Time step number
 * @param finest_level Finest AMR level
 * @param varnames Vector of variable names to output
 * @param mf Vector of MultiFabs containing the data (already gathered)
 * @param geom Vector of geometry objects for each level
 * @param time Current simulation time
 * @param istep Vector of time step indices for each level
 * @param ref_ratio Refinement ratios between levels
 */
void writePlotfile(const std::string& plot_file,
                   const int step,
                   const int finest_level,
                   const amrex::Vector<std::string>& varnames,
                   const amrex::Vector<amrex::MultiFab>& mf,
                   const amrex::Vector<amrex::Geometry>& geom,
                   const amrex::Real time,
                   const amrex::Vector<int>& istep,
                   const amrex::Vector<amrex::IntVect>& ref_ratio);

/**
 * @brief Write Tecplot ASCII format file with surface data (eta, u, v) - single level
 *
 * @param plot_file Base name for the output file
 * @param step Current time step number
 * @param geom Geometry object
 * @param eta Surface elevation MultiFab
 * @param u Horizontal velocity (x-direction) MultiFab
 * @param v Horizontal velocity (y-direction) MultiFab
 * @param ba_root BoxArray for root processor
 * @param dmap_root DistributionMapping for root processor
 * @param time Current simulation time
 */
void writeTecplot(const std::string& plot_file,
                  const int step,
                  const amrex::Geometry& geom,
                  const amrex::MultiFab& eta,
                  const amrex::MultiFab& u,
                  const amrex::MultiFab& v,
                  const amrex::BoxArray& ba_root,
                  const amrex::DistributionMapping& dmap_root,
                  const amrex::Real time);

/**
 * @brief Write multi-level Tecplot ASCII format file with surface data for nested AMR
 *
 * Writes all AMR levels to a single Tecplot file with multiple zones.
 * Each level is written as a separate zone with appropriate grid resolution.
 *
 * @param plot_file Base name for the output file
 * @param step Current time step number
 * @param finest_level Finest AMR level to output
 * @param geom Vector of Geometry objects for each level
 * @param eta Vector of surface elevation MultiFabs
 * @param u Vector of horizontal velocity (x) MultiFabs
 * @param v Vector of horizontal velocity (y) MultiFabs
 * @param ba_root Vector of BoxArray for root processor
 * @param dmap_root Vector of DistributionMapping for root processor
 * @param time Current simulation time
 */
void writeTecplotMultiLevel(const std::string& plot_file,
                             const int step,
                             const int finest_level,
                             const amrex::Vector<amrex::Geometry>& geom,
                             const amrex::Vector<amrex::MultiFab>& eta,
                             const amrex::Vector<amrex::MultiFab>& u,
                             const amrex::Vector<amrex::MultiFab>& v,
                             const amrex::Vector<amrex::BoxArray>& ba_root,
                             const amrex::Vector<amrex::DistributionMapping>& dmap_root,
                             const amrex::Real time);

/**
 * @brief Write multi-level Tecplot file with CRESCENT diagnostic variables
 *
 * Extended version that includes CRESCENT diagnostic variables for nested AMR simulations.
 * CRESCENT variables are only output for refined levels (lev >= 1).
 *
 * @param plot_file Base name for the output file
 * @param step Current time step number
 * @param finest_level Finest AMR level to output
 * @param geom Vector of Geometry objects for each level
 * @param eta Vector of surface elevation MultiFabs
 * @param u Vector of horizontal velocity (x) MultiFabs
 * @param v Vector of horizontal velocity (y) MultiFabs
 * @param ba_root Vector of BoxArray for root processor
 * @param dmap_root Vector of DistributionMapping for root processor
 * @param time Current simulation time
 * @param diagnostic_levels Grouped CRESCENT output sources:
 *        - `core_surface_derivatives` are written on level 0
 *        - `nested_output_derivatives` are written on refined levels
 *        - `eta_x_core_reference_l1` / `phi_x_core_reference_l1` remain the
 *          L0-on-L1 comparison fields for refined-level diagnostics
 */
void writeTecplotMultiLevelCRESCENT(const std::string& plot_file,
                                   const int step,
                                   const int finest_level,
                                   const amrex::Vector<amrex::Geometry>& geom,
                                   const amrex::Vector<amrex::MultiFab>& eta,
                                   const amrex::Vector<amrex::MultiFab>& phi,
                                   const amrex::Vector<amrex::MultiFab>& u,
                                   const amrex::Vector<amrex::MultiFab>& v,
                                   const amrex::Vector<amrex::BoxArray>& ba_root,
                                   const amrex::Vector<amrex::DistributionMapping>& dmap_root,
                                   const amrex::Real time,
                                   const NestedTecplotDiagnosticLevels& diagnostic_levels);

/**
 * @brief Write multi-level Tecplot movie file (eta only, L1 restricted to Ω_A)
 *
 * Optimized for animation/visualization. Outputs only surface elevation (η) with:
 * - L0: Full domain
 * - L1: Only the Ω_A (accurate) region where CRESCENT provides highest quality results
 *
 * This produces smaller files focused on the high-fidelity region, ideal for movies.
 *
 * @param plot_file Base name for the output file
 * @param step Current time step number
 * @param finest_level Finest AMR level to output
 * @param geom Vector of Geometry objects for each level
 * @param eta Vector of surface elevation MultiFabs
 * @param ba_root Vector of BoxArray for root processor
 * @param dmap_root Vector of DistributionMapping for root processor
 * @param time Current simulation time
 * @param crescent_domains Vector of CRESCENT domain structures (contains box_A for Ω_A)
 */
void writeTecplotMovie(const std::string& plot_file,
                       const int step,
                       const int finest_level,
                       const amrex::Vector<amrex::Geometry>& geom,
                       const amrex::Vector<amrex::MultiFab>& eta,
                       const amrex::Vector<amrex::BoxArray>& ba_root,
                       const amrex::Vector<amrex::DistributionMapping>& dmap_root,
                       const amrex::Real time,
                       const amrex::Vector<crescent::nested::domains::DomainStructure>& crescent_domains);

} // namespace crescent::io

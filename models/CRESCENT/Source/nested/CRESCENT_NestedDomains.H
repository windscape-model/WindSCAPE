#ifndef CRESCENT_NESTED_DOMAINS_H
#define CRESCENT_NESTED_DOMAINS_H

#include <AMReX_Geometry.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_TagBox.H>
#include <AMReX_Vector.H>

#include "CRESCENT_MeshHierarchyView.H"

namespace crescent::nested {

struct AlgorithmConfig;
struct RuntimeState;

}

/**
 * @file CRESCENT_NestedDomains.H
 * @brief Domain structure definitions for the nested CRESCENT path
 *
 * CRESCENT uses nested subdomains within the L1 grid:
 *   Ω_A (Accurate):  Pure L1 nested dynamics, centered region
 *   Ω_S (Spectral):  Region where FFT results are reliable
 *   Ω_M (Mixing):    Ω_S \ Ω_A, smooth nested↔core transition
 *   Ω_E (Exterior):  L1 \ Ω_S, pure core interpolation
 */

namespace crescent::nested::domains {

struct CouplingZoneSetupContext
{
    const crescent::nested::AlgorithmConfig& config;
    crescent::nested::RuntimeState& runtime;
    int verbosity;
    crescent::shared::MeshHierarchyView mesh;

    crescent::shared::LevelMeshView levelMesh(int lev) const
    {
        return mesh.levelMesh(lev);
    }
};

void initializeLevelCouplingZones(
    int lev,
    const CouplingZoneSetupContext& context);

struct L1RefinementRegion
{
    amrex::Real x_lo;
    amrex::Real x_hi;
    amrex::Real y_lo;
    amrex::Real y_hi;
    amrex::Real dx_x;
    amrex::Real dx_y;
    amrex::Real Lx;
    amrex::Real Ly;
    int nx;
    int ny;
};

struct L1RefinementTaggingContext
{
    const crescent::nested::AlgorithmConfig& config;
    const amrex::Geometry& geom;
    int verbosity;
    amrex::IntVect ref_ratio;
    amrex::IntVect blocking_factor;
};

void tagL1RefinementRegion(
    amrex::TagBoxArray& tags,
    const L1RefinementTaggingContext& context);

/**
 * @struct DomainStructure
 * @brief Stores boundaries for CRESCENT subdomains
 *
 * All boundaries are in physical coordinates and index space.
 */
struct DomainStructure {
    // Configuration parameters
    amrex::Real accurate_fraction;   ///< Margin for Ω_A (default: 0.25 = 25%)
    amrex::Real spectral_fraction;   ///< Margin for Ω_S (default: 0.15 = 15%)

    // Physical boundaries - L1 refined box extent
    amrex::Real x_lo_L1, x_hi_L1;
    amrex::Real y_lo_L1, y_hi_L1;

    // Physical boundaries - Accurate domain Ω_A
    amrex::Real x_lo_A, x_hi_A;
    amrex::Real y_lo_A, y_hi_A;

    // Physical boundaries - Spectral subdomain Ω_S
    amrex::Real x_lo_S, x_hi_S;
    amrex::Real y_lo_S, y_hi_S;

    // Index space boundaries
    amrex::Box box_A;                ///< Index space box for Ω_A
    amrex::Box box_S;                ///< Index space box for Ω_S

    // Zone widths (for mask generation)
    amrex::Real W_exterior_x, W_exterior_y;  ///< Width of exterior zones (Ω_E) in x,y
    amrex::Real W_mix_x, W_mix_y;            ///< Width of mixing zones (Ω_M) in x,y

    // Domain info
    int ni, nj;                      ///< Grid dimensions
    amrex::Real Lx, Ly;              ///< Physical domain size
};

/**
 * @brief Compute domain structure for CRESCENT algorithm
 *
 * Given L1 geometry, grid extent, and margin fractions, computes physical and
 * index-space boundaries for all CRESCENT subdomains.
 *
 * Ω_S and Ω_A are defined relative to the L1 refined box, not the problem
 * domain, so every zone lies inside L1.
 *
 * @param geom              L1 level geometry
 * @param grids             L1 BoxArray (to get actual refined box extent)
 * @param accurate_fraction Margin for Ω_A as fraction of L1 box (0.25 = 25% margin)
 * @param spectral_fraction Margin for Ω_S as fraction of L1 box (0.15 = 15% margin)
 * @return DomainStructure with all boundaries computed
 *
 * Example: For L1 box [40, 160]m with accurate_fraction=0.25:
 *   - L1 extent: Lx=120m
 *   - Ω_A: [40+30, 160-30] = [70, 130] (centered 50% of L1 box)
 *   - Margin: 30m on each side (25% of 120m)
 */
DomainStructure computeDomainStructure(
    const amrex::Geometry& geom,
    const amrex::BoxArray& grids,
    amrex::Real accurate_fraction = 0.25,
    amrex::Real spectral_fraction = 0.15
);

/**
 * @brief Print domain-structure diagnostics
 */
void printDomainStructure(const DomainStructure& domain);

} // namespace crescent::nested::domains

#endif // CRESCENT_NESTED_DOMAINS_H

#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>


namespace crescent::interpolation {

/**
 * @brief Spatial interpolation from parent (L0) to child (L1) grid
 *
 * Uses custom 2D interpolation (bilinear or quadratic) for surface state fields.
 * This is used to provide initial conditions for the fine level from the coarse
 * parent solution.
 *
 * Key features:
 * - Spatial interpolation ONLY (single time level)
 * - Uses single-box coarse patch pattern (MPI-safe, multi-box safe)
 * - nz=1 compatible (2D physics in 3D build)
 * - Strategy: ParallelCopy gathers coarse data → interpolate → FillBoundary
 *
 * @param parent_eta Coarse grid surface elevation η
 * @param parent_phi Coarse grid velocity potential ψ
 * @param child_eta Fine grid surface elevation (output)
 * @param child_phi Fine grid velocity potential (output)
 * @param parent_geom Coarse grid geometry
 * @param child_geom Fine grid geometry
 * @param ref_ratio Refinement ratio between levels
 * @param coarse_patch_scratch Caller-owned coarse-grid patch scratch
 * @param interpolation_order Interpolation order (1=bilinear, 2=quadratic, default=2)
 */
void interpolateParentToChild(
    const amrex::MultiFab& parent_eta,
    const amrex::MultiFab& parent_phi,
    amrex::MultiFab& child_eta,
    amrex::MultiFab& child_phi,
    const amrex::Geometry& parent_geom,
    const amrex::Geometry& child_geom,
    const amrex::IntVect& ref_ratio,
    amrex::MultiFab& coarse_patch_scratch,
    int interpolation_order = 2
);

/**
 * @brief Spatial interpolation from one parent field to one child field
 *
 * Uses the same custom 2D interpolation machinery as interpolateParentToChild(),
 * but for scalar fields such as CRESCENT derivative components.
 *
 * @param parent Coarse-grid scalar field
 * @param child Fine-grid scalar field (output)
 * @param parent_geom Coarse-grid geometry
 * @param child_geom Fine-grid geometry
 * @param ref_ratio Refinement ratio between levels
 * @param coarse_patch_scratch Caller-owned coarse-grid patch scratch
 * @param interpolation_order Interpolation order (1=bilinear, 2=quadratic, default=2)
 */
void interpolateScalarParentToChild(
    const amrex::MultiFab& parent,
    amrex::MultiFab& child,
    const amrex::Geometry& parent_geom,
    const amrex::Geometry& child_geom,
    const amrex::IntVect& ref_ratio,
    amrex::MultiFab& coarse_patch_scratch,
    int interpolation_order = 2
);

} // namespace crescent::interpolation

#ifndef CRESCENT_NESTED_C3PERIODICITY_H_
#define CRESCENT_NESTED_C3PERIODICITY_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include "CRESCENT_NestedDomains.H"
#include "CRESCENT_NestedTypes.H"

/**
 * @file CRESCENT_NestedC3Periodicity.H
 * @brief C3 interface repair of the PTF periodic copies.
 *
 * The repair adds a correction sum_r alpha_r Q_r(xi), Q_r = (1-xi)^4 xi^r,
 * inside the buffer next to each child box edge, so that the value and the
 * first three derivatives match across the periodic wrap. Omega_A and the
 * mixing annulus are not modified. The coefficients come from a regularized
 * least-squares fit (boundary penalty lambda_b, coefficient regularization
 * lambda_c); the basis is computed once per level (computeC3InterfaceBasis).
 */

namespace crescent::nested::c3 {

inline constexpr int InterfaceBasisOrder = 4;

struct C3InterfaceRepairOptions {
    int buffer_width_x = 12;
    int buffer_width_y = 12;
    amrex::Real boundary_penalty_weight = 0.01;
    amrex::Real coefficient_regularization_weight = 10.0;
};

/**
 * @brief Precompute the C3 interface basis (Qx, Qy, Bx, By) into a workspace.
 *
 * The basis depends only on the deterministic grid parameters
 * (InterfaceBasisOrder, nx, ny, buffer widths bx/by, and cell sizes dx/dy).
 * Call once at level setup; applyC3InterfaceSeamRepair then reads the cached
 * basis and does no per-call basis construction.
 */
void computeC3InterfaceBasis(
    crescent::nested::C3RepairWorkspace& workspace,
    int k,
    int nx,
    int ny,
    int bx,
    int by,
    amrex::Real dx,
    amrex::Real dy);

/**
 * @brief Repair one field: an x pass on every row, then a y pass on every
 * column. Each line gets the coefficients of a regularized least-squares fit
 * of the four wrap jumps (value and first three derivatives).
 */
void applyC3InterfaceSeamRepair(
    amrex::MultiFab& field,
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    const C3InterfaceRepairOptions& options,
    crescent::nested::C3RepairWorkspace& workspace);

} // namespace crescent::nested::c3

#endif // CRESCENT_NESTED_C3PERIODICITY_H_

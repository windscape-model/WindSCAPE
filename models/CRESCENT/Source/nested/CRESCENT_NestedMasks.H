#ifndef CRESCENT_NESTED_MASKS_H
#define CRESCENT_NESTED_MASKS_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_NestedDomains.H"

namespace crescent::nested::masks {

/**
 * @brief Build the blending function sigma of the child level.
 *
 * sigma = 1 - (1 - S(d_x / W_x)) (1 - S(d_y / W_y)), where S is the septic
 * smoothstep, d_x and d_y are the distances outside Omega_A, and W_x, W_y
 * are the widths of the mixing annulus between Omega_A and Omega_S. Hence
 * sigma = 0 in Omega_A, sigma = 1 on and outside the Omega_S boundary, and
 * sigma is C^3 everywhere, the corner seams included.
 */
void generateBlendingFunction(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain);

} // namespace crescent::nested::masks

#endif // CRESCENT_NESTED_MASKS_H

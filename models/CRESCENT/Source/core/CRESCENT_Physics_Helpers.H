#pragma once

#include <AMReX_MultiFab.H>

/**
 * @file CRESCENT_Physics_Helpers.H
 * @brief Production helpers for the active core DNO path
 *
 * Keep this surface limited to helpers that are actually used by the live
 * solver. The broader experimental DNO helper decomposition was retired from
 * active Source once `crescent::dno::detail` became the real shared helper layer.
 */

namespace crescent::physics::helpers {

/**
 * @brief Compute η^k/k! for k=0,1,...,M-1
 *
 * Used in Taylor expansion: φ(x,y,η) = Σ η^k/k! * ∂^k φ/∂z^k |_{z=0}
 *
 * @param eta_pow Output: η^k/k! for each order
 * @param eta Input: surface elevation
 * @param perturbation_order M (number of terms in expansion)
 */
void computeEtaPowersWithFactorial(
    amrex::MultiFab& eta_pow,
    const amrex::MultiFab& eta,
    int perturbation_order);

} // namespace crescent::physics::helpers

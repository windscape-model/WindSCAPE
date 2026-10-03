#pragma once

#include <AMReX_REAL.H>

/**
 * @file CRESCENT_DNOConfig.H
 * @brief Live Dirichlet-to-Neumann operator configuration
 */

namespace crescent::config {

/**
 * @brief Live configuration for Dirichlet-to-Neumann Operator (DNO)
 *
 * Keep only fields that are read on the active DNO path.
 */
struct DNOConfig {
    int perturbation_order = 3;
    int verbosity = 1;
};

inline DNOConfig makeDNOConfig(int perturbation_order, int verbosity)
{
    return {perturbation_order, verbosity};
}

} // namespace crescent::config

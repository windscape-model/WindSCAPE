#ifndef CRESCENT_NESTED_C3_LINEAR_SOLVERS_H
#define CRESCENT_NESTED_C3_LINEAR_SOLVERS_H

#include <array>
#include <vector>

#include <AMReX_Real.H>

namespace crescent::nested::c3 {

bool solve_least_squares(
    const std::vector<std::vector<amrex::Real>>& B,
    const amrex::Real Delta[4],
    int k,
    amrex::Real alpha[],
    amrex::Real boundary_penalty_weight = 0.0,
    amrex::Real target_left = 0.0,
    amrex::Real target_right = 0.0,
    amrex::Real coefficient_regularization = 1.0);

bool solve_coupled_line_system(
    const std::vector<std::vector<amrex::Real>>& B,
    const std::vector<std::array<amrex::Real, 4>>& deltas,
    const std::vector<amrex::Real>& target_left,
    const std::vector<amrex::Real>& target_right,
    int ncoef,
    amrex::Real boundary_penalty_weight,
    amrex::Real coefficient_regularization,
    amrex::Real transverse_regularization_weight,
    std::vector<amrex::Real>& alpha_flat);

bool solve_mean_coupled_line_system(
    const std::vector<std::vector<amrex::Real>>& B,
    const std::vector<std::array<amrex::Real, 4>>& deltas,
    const std::vector<amrex::Real>& target_left,
    const std::vector<amrex::Real>& target_right,
    int ncoef,
    amrex::Real boundary_penalty_weight,
    amrex::Real coefficient_regularization,
    amrex::Real mean_regularization_weight,
    std::vector<amrex::Real>& alpha_flat);

void build_periodic_real_basis(
    int npts,
    int nharmonics,
    std::vector<std::vector<amrex::Real>>& basis);

void project_coefficients_onto_basis(
    std::vector<amrex::Real>& coeff_flat,
    int nlines,
    int ncoef,
    const std::vector<std::vector<amrex::Real>>& basis);

bool solve_modal_edge_system(
    const std::vector<std::vector<amrex::Real>>& B,
    const std::vector<std::array<amrex::Real, 4>>& deltas,
    const std::vector<amrex::Real>& target_left,
    const std::vector<amrex::Real>& target_right,
    const std::vector<std::vector<amrex::Real>>& transverse_basis,
    int ncoef,
    amrex::Real boundary_penalty_weight,
    amrex::Real coefficient_regularization,
    std::vector<amrex::Real>& alpha_flat);

} // namespace crescent::nested::c3

#endif // CRESCENT_NESTED_C3_LINEAR_SOLVERS_H

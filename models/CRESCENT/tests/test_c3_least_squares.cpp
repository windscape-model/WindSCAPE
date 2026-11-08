#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "CRESCENT_NestedC3LinearSolvers.H"

namespace {

// The C3 repair fits 2k = 8 coefficients per line to the four wrap jumps.
// Check that solve_least_squares returns the solution of the regularized,
// endpoint-penalized normal system that it is documented to solve.
int runTest()
{
    using amrex::Real;

    constexpr int ncoef = 8;
    std::vector<std::vector<Real>> B(4, std::vector<Real>(ncoef));
    for (int n = 0; n < 4; ++n) {
        for (int i = 0; i < ncoef; ++i) {
            B[n][i] = std::sin(1.0 + 0.7 * n + 1.3 * i) * std::pow(2.0, n);
        }
    }
    const Real Delta[4] = {0.4, -0.2, 0.15, 0.3};
    const Real target_left = 0.05;
    const Real target_right = -0.08;
    const Real boundary_penalty = 0.01;
    const Real coefficient_regularization = 10.0;

    Real alpha[ncoef];
    if (!crescent::nested::c3::solve_least_squares(
            B, Delta, ncoef, alpha, boundary_penalty, target_left,
            target_right, coefficient_regularization)) {
        std::cerr << "C3 least-squares solve reported a singular system\n";
        return 1;
    }

    // Normal system: (B^T B + lambda I + mu^2 E) alpha = -B^T Delta + mu^2 t,
    // with lambda = (1e-8 + lambda_c) s, mu^2 = lambda_b s, s = tr(B^T B)/k.
    std::vector<std::vector<Real>> M(ncoef, std::vector<Real>(ncoef, 0.0));
    std::vector<Real> rhs(ncoef, 0.0);
    Real scale = 0.0;
    for (int i = 0; i < ncoef; ++i) {
        for (int j = 0; j < ncoef; ++j) {
            for (int n = 0; n < 4; ++n) {
                M[i][j] += B[n][i] * B[n][j];
            }
        }
        for (int n = 0; n < 4; ++n) {
            rhs[i] -= B[n][i] * Delta[n];
        }
        scale += M[i][i];
    }
    scale = scale / ncoef + 1.0e-14;
    for (int i = 0; i < ncoef; ++i) {
        M[i][i] += (1.0e-8 + coefficient_regularization) * scale;
    }
    const Real mu_sq = boundary_penalty * scale;
    M[0][0] += mu_sq;
    rhs[0] += mu_sq * target_left;
    M[ncoef / 2][ncoef / 2] += mu_sq;
    rhs[ncoef / 2] += mu_sq * target_right;

    Real residual = 0.0;
    Real rhs_norm = 0.0;
    for (int i = 0; i < ncoef; ++i) {
        Real row = -rhs[i];
        for (int j = 0; j < ncoef; ++j) {
            row += M[i][j] * alpha[j];
        }
        residual = std::max(residual, std::abs(row));
        rhs_norm = std::max(rhs_norm, std::abs(rhs[i]));
    }
    if (residual > 1.0e-12 * rhs_norm) {
        std::cerr << "C3 least-squares relative residual = "
                  << residual / rhs_norm << '\n';
        return 1;
    }
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = amrex::ParallelDescriptor::IOProcessor() ? runTest() : 0;
    amrex::Finalize();
    return result;
}

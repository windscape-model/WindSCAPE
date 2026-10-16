#ifndef CRESCENT_NESTED_C3_LINEAR_SOLVERS_H
#define CRESCENT_NESTED_C3_LINEAR_SOLVERS_H

#include <vector>

#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>

namespace crescent::nested::c3 {

// Regularized least-squares fit of the C3 repair coefficients on fixed-size
// arrays, callable from the batched GPU kernel and (through
// solve_least_squares below) from the host path, so both paths run the same
// arithmetic. B is the 4 x k basis-jump matrix (k <= 8), Delta the four seam
// jumps; alpha receives the k coefficients. Returns false if the normal
// system is singular.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool solve_least_squares_device(
    const amrex::Real B[4][8],
    const amrex::Real Delta[4],
    int k,
    amrex::Real alpha[],
    amrex::Real boundary_penalty_weight,
    amrex::Real target_left,
    amrex::Real target_right,
    amrex::Real coefficient_regularization)
{
    using amrex::Real;
    Real M[8][8] = {};
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            Real sum = 0.0;
            for (int n = 0; n < 4; ++n) {
                sum += B[n][i] * B[n][j];
            }
            M[i][j] = sum;
        }
    }

    Real scale = 0.0;
    for (int i = 0; i < k; ++i) {
        scale += M[i][i];
    }
    scale = (scale / Real(k)) + 1.0e-14;

    const Real lambda_sq = 1.0e-8 * scale + coefficient_regularization * scale;
    for (int i = 0; i < k; ++i) {
        M[i][i] += lambda_sq;
    }

    Real rhs[8] = {};
    for (int i = 0; i < k; ++i) {
        Real sum = 0.0;
        for (int n = 0; n < 4; ++n) {
            sum += B[n][i] * Delta[n];
        }
        rhs[i] = -sum;
    }

    // Endpoint penalty (see the host build_least_squares_system): +mu^2 t.
    if (boundary_penalty_weight > 0.0) {
        const Real mu_sq = boundary_penalty_weight * scale;
        M[0][0] += mu_sq;
        rhs[0] += mu_sq * target_left;
        const int idx_right = k / 2;
        M[idx_right][idx_right] += mu_sq;
        rhs[idx_right] += mu_sq * target_right;
    }

    Real aug[8][9] = {};
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            aug[i][j] = M[i][j];
        }
        aug[i][k] = rhs[i];
    }

    for (int p = 0; p < k; ++p) {
        int pivot = p;
        Real max_val = amrex::Math::abs(aug[p][p]);
        for (int i = p + 1; i < k; ++i) {
            const Real a = amrex::Math::abs(aug[i][p]);
            if (a > max_val) {
                max_val = a;
                pivot = i;
            }
        }
        if (max_val < 1.0e-14) {
            return false;
        }
        if (pivot != p) {
            for (int j = 0; j < k + 1; ++j) {
                const Real tmp = aug[p][j];
                aug[p][j] = aug[pivot][j];
                aug[pivot][j] = tmp;
            }
        }
        for (int i = p + 1; i < k; ++i) {
            const Real factor = aug[i][p] / aug[p][p];
            for (int j = p; j < k + 1; ++j) {
                aug[i][j] -= factor * aug[p][j];
            }
        }
    }

    for (int i = k - 1; i >= 0; --i) {
        alpha[i] = aug[i][k];
        for (int j = i + 1; j < k; ++j) {
            alpha[i] -= aug[i][j] * alpha[j];
        }
        alpha[i] /= aug[i][i];
    }
    return true;
}

// Host entry point for basis matrices held as nested vectors (the CPU repair
// path). It copies B into the fixed-size layout and runs the same routine.
inline bool solve_least_squares(
    const std::vector<std::vector<amrex::Real>>& B,
    const amrex::Real Delta[4],
    int k,
    amrex::Real alpha[],
    amrex::Real boundary_penalty_weight,
    amrex::Real target_left,
    amrex::Real target_right,
    amrex::Real coefficient_regularization)
{
    AMREX_ALWAYS_ASSERT(k > 0 && k <= 8);
    amrex::Real Bfix[4][8] = {};
    for (int n = 0; n < 4; ++n) {
        for (int i = 0; i < k; ++i) {
            Bfix[n][i] = B[n][i];
        }
    }
    return solve_least_squares_device(
        Bfix, Delta, k, alpha, boundary_penalty_weight,
        target_left, target_right, coefficient_regularization);
}

} // namespace crescent::nested::c3

#endif // CRESCENT_NESTED_C3_LINEAR_SOLVERS_H

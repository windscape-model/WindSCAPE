#include "CRESCENT_NestedC3LinearSolvers.H"

// Shared linear-system and modal solver kernels for nested C3 repair.

#include <AMReX_BLassert.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

namespace crescent::nested::c3 {

using namespace amrex;

namespace {

void build_least_squares_system(
    const std::vector<std::vector<Real>>& B,
    const Real Delta[4],
    int k,
    std::vector<std::vector<Real>>& M,
    std::vector<Real>& rhs,
    Real boundary_penalty_weight,
    Real target_left,
    Real target_right,
    Real coefficient_regularization,
    Real& scale_out)
{
    M.assign(k, std::vector<Real>(k, 0.0));
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
    scale_out = scale;

    const Real base_lambda_sq = 1.0e-8 * scale;
    const Real coef_lambda_sq = coefficient_regularization * scale;
    const Real lambda_sq = base_lambda_sq + coef_lambda_sq;
    for (int i = 0; i < k; ++i) {
        M[i][i] += lambda_sq;
    }

    rhs.assign(k, 0.0);
    for (int i = 0; i < k; ++i) {
        Real sum = 0.0;
        for (int n = 0; n < 4; ++n) {
            sum += B[n][i] * Delta[n];
        }
        rhs[i] = -sum;
    }

    if (boundary_penalty_weight > 0.0) {
        const Real mu_sq = boundary_penalty_weight * scale;
        M[0][0] += mu_sq;
        rhs[0] -= mu_sq * target_left;

        const int idx_right = k / 2;
        M[idx_right][idx_right] += mu_sq;
        rhs[idx_right] -= mu_sq * target_right;
    }
}

bool solve_dense_system(
    std::vector<std::vector<Real>>& A,
    std::vector<Real>& b)
{
    const int n = static_cast<int>(A.size());
    AMREX_ALWAYS_ASSERT(n == static_cast<int>(b.size()));

    for (int p = 0; p < n; ++p) {
        int pivot = p;
        Real max_val = std::abs(A[p][p]);
        for (int i = p + 1; i < n; ++i) {
            if (std::abs(A[i][p]) > max_val) {
                max_val = std::abs(A[i][p]);
                pivot = i;
            }
        }

        if (max_val < 1.0e-14) {
            return false;
        }

        if (pivot != p) {
            std::swap(A[p], A[pivot]);
            std::swap(b[p], b[pivot]);
        }

        for (int i = p + 1; i < n; ++i) {
            const Real factor = A[i][p] / A[p][p];
            for (int j = p; j < n; ++j) {
                A[i][j] -= factor * A[p][j];
            }
            b[i] -= factor * b[p];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            b[i] -= A[i][j] * b[j];
        }
        b[i] /= A[i][i];
    }

    return true;
}

bool solve_dense_multiple_rhs(
    const std::vector<std::vector<Real>>& A_in,
    const std::vector<std::vector<Real>>& B_in,
    std::vector<std::vector<Real>>& X_out)
{
    const int n = static_cast<int>(A_in.size());
    AMREX_ALWAYS_ASSERT(n == static_cast<int>(B_in.size()));
    X_out.assign(n, std::vector<Real>(n, 0.0));

    for (int col = 0; col < n; ++col) {
        auto A = A_in;
        std::vector<Real> x(n, 0.0);
        for (int row = 0; row < n; ++row) {
            x[row] = B_in[row][col];
        }
        if (!solve_dense_system(A, x)) {
            return false;
        }
        for (int row = 0; row < n; ++row) {
            X_out[row][col] = x[row];
        }
    }

    return true;
}

} // namespace

bool solve_least_squares(
    const std::vector<std::vector<Real>>& B,
    const Real Delta[4],
    int k,
    Real alpha[],
    Real boundary_penalty_weight,
    Real target_left,
    Real target_right,
    Real coefficient_regularization)
{
    constexpr int max_line_coefficients = 8;
    AMREX_ALWAYS_ASSERT(k > 0);
    AMREX_ALWAYS_ASSERT(k <= max_line_coefficients);

    std::array<std::array<Real, max_line_coefficients>, max_line_coefficients> M{};
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

    const Real base_lambda_sq = 1.0e-8 * scale;
    const Real coef_lambda_sq = coefficient_regularization * scale;
    const Real lambda_sq = base_lambda_sq + coef_lambda_sq;
    for (int i = 0; i < k; ++i) {
        M[i][i] += lambda_sq;
    }

    std::array<Real, max_line_coefficients> rhs{};
    for (int i = 0; i < k; ++i) {
        Real sum = 0.0;
        for (int n = 0; n < 4; ++n) {
            sum += B[n][i] * Delta[n];
        }
        rhs[i] = -sum;
    }

    if (boundary_penalty_weight > 0.0) {
        const Real mu_sq = boundary_penalty_weight * scale;

        M[0][0] += mu_sq;
        rhs[0] -= mu_sq * target_left;

        const int idx_right = k / 2;
        M[idx_right][idx_right] += mu_sq;
        rhs[idx_right] -= mu_sq * target_right;
    }

    std::array<std::array<Real, max_line_coefficients + 1>, max_line_coefficients> aug{};
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            aug[i][j] = M[i][j];
        }
        aug[i][k] = rhs[i];
    }

    for (int p = 0; p < k; ++p) {
        int pivot = p;
        Real max_val = std::abs(aug[p][p]);
        for (int i = p + 1; i < k; ++i) {
            if (std::abs(aug[i][p]) > max_val) {
                max_val = std::abs(aug[i][p]);
                pivot = i;
            }
        }

        if (max_val < 1.0e-14) {
            return false;
        }

        if (pivot != p) {
            std::swap(aug[p], aug[pivot]);
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

    Real alpha_max = 0.0;
    for (int i = 0; i < k; ++i) {
        alpha_max = std::max(alpha_max, std::abs(alpha[i]));
    }

    if (alpha_max > 50.0) {
        // Large coefficients detected; the caller may want stronger
        // coefficient regularization for this case.
    }

    return true;
}

bool solve_coupled_line_system(
    const std::vector<std::vector<Real>>& B,
    const std::vector<std::array<Real, 4>>& deltas,
    const std::vector<Real>& target_left,
    const std::vector<Real>& target_right,
    int ncoef,
    Real boundary_penalty_weight,
    Real coefficient_regularization,
    Real transverse_regularization_weight,
    std::vector<Real>& alpha_flat)
{
    const int nlines = static_cast<int>(deltas.size());
    AMREX_ALWAYS_ASSERT(nlines == static_cast<int>(target_left.size()));
    AMREX_ALWAYS_ASSERT(nlines == static_cast<int>(target_right.size()));

    std::vector<std::vector<std::vector<Real>>> diag(nlines);
    std::vector<std::vector<Real>> rhs(nlines);

    Real scale = 0.0;
    for (int line = 0; line < nlines; ++line) {
        Real line_delta[4] = {
            deltas[line][0], deltas[line][1], deltas[line][2], deltas[line][3]
        };
        Real line_scale = 0.0;
        build_least_squares_system(
            B, line_delta, ncoef, diag[line], rhs[line],
            boundary_penalty_weight, target_left[line], target_right[line],
            coefficient_regularization, line_scale);
        scale = std::max(scale, line_scale);
    }

    if (nlines == 1 || transverse_regularization_weight <= 0.0) {
        alpha_flat.assign(nlines * ncoef, 0.0);
        for (int line = 0; line < nlines; ++line) {
            auto A = diag[line];
            auto x = rhs[line];
            if (!solve_dense_system(A, x)) {
                return false;
            }
            for (int r = 0; r < ncoef; ++r) {
                alpha_flat[line * ncoef + r] = x[r];
            }
        }
        return true;
    }

    const Real nu = transverse_regularization_weight * scale;
    std::vector<std::vector<Real>> offdiag_template(
        ncoef, std::vector<Real>(ncoef, 0.0));
    for (int r = 0; r < ncoef; ++r) {
        offdiag_template[r][r] = -nu;
    }

    for (int line = 0; line < nlines; ++line) {
        const Real diag_add = (line == 0 || line == nlines - 1) ? nu : (2.0 * nu);
        for (int r = 0; r < ncoef; ++r) {
            diag[line][r][r] += diag_add;
        }
    }

    std::vector<std::vector<std::vector<Real>>> upper(
        nlines - 1, offdiag_template);

    for (int line = 1; line < nlines; ++line) {
        std::vector<std::vector<Real>> solved_upper;
        if (!solve_dense_multiple_rhs(diag[line - 1], upper[line - 1], solved_upper)) {
            return false;
        }

        auto solved_rhs = rhs[line - 1];
        auto diag_copy = diag[line - 1];
        if (!solve_dense_system(diag_copy, solved_rhs)) {
            return false;
        }

        for (int i = 0; i < ncoef; ++i) {
            for (int j = 0; j < ncoef; ++j) {
                diag[line][i][j] -= (-nu) * solved_upper[i][j];
            }
            rhs[line][i] -= (-nu) * solved_rhs[i];
        }

        upper[line - 1] = std::move(solved_upper);
    }

    std::vector<std::vector<Real>> solution(
        nlines, std::vector<Real>(ncoef, 0.0));
    {
        auto last_diag = diag[nlines - 1];
        auto last_rhs = rhs[nlines - 1];
        if (!solve_dense_system(last_diag, last_rhs)) {
            return false;
        }
        solution[nlines - 1] = std::move(last_rhs);
    }

    for (int line = nlines - 2; line >= 0; --line) {
        auto rhs_eff = rhs[line];
        for (int i = 0; i < ncoef; ++i) {
            for (int j = 0; j < ncoef; ++j) {
                rhs_eff[i] -= upper[line][i][j] * solution[line + 1][j];
            }
        }
        auto diag_copy = diag[line];
        if (!solve_dense_system(diag_copy, rhs_eff)) {
            return false;
        }
        solution[line] = std::move(rhs_eff);
    }

    alpha_flat.assign(nlines * ncoef, 0.0);
    for (int line = 0; line < nlines; ++line) {
        for (int r = 0; r < ncoef; ++r) {
            alpha_flat[line * ncoef + r] = solution[line][r];
        }
    }

    return true;
}

bool solve_mean_coupled_line_system(
    const std::vector<std::vector<Real>>& B,
    const std::vector<std::array<Real, 4>>& deltas,
    const std::vector<Real>& target_left,
    const std::vector<Real>& target_right,
    int ncoef,
    Real boundary_penalty_weight,
    Real coefficient_regularization,
    Real mean_regularization_weight,
    std::vector<Real>& alpha_flat)
{
    const int nlines = static_cast<int>(deltas.size());
    AMREX_ALWAYS_ASSERT(nlines == static_cast<int>(target_left.size()));
    AMREX_ALWAYS_ASSERT(nlines == static_cast<int>(target_right.size()));

    alpha_flat.assign(nlines * ncoef, 0.0);
    if (nlines == 0) {
        return true;
    }

    Real zero_delta[4] = {0.0, 0.0, 0.0, 0.0};
    std::vector<std::vector<Real>> common_A;
    std::vector<Real> rhs_dummy;
    Real scale = 0.0;
    build_least_squares_system(
        B, zero_delta, ncoef, common_A, rhs_dummy,
        boundary_penalty_weight, 0.0, 0.0,
        coefficient_regularization, scale);

    std::vector<std::vector<Real>> rhs_lines(
        nlines, std::vector<Real>(ncoef, 0.0));
    std::vector<Real> mean_rhs(ncoef, 0.0);
    for (int line = 0; line < nlines; ++line) {
        for (int i = 0; i < ncoef; ++i) {
            Real sum = 0.0;
            for (int n = 0; n < 4; ++n) {
                sum += B[n][i] * deltas[line][n];
            }
            rhs_lines[line][i] = -sum;
        }

        if (boundary_penalty_weight > 0.0) {
            const Real mu_sq = boundary_penalty_weight * scale;
            rhs_lines[line][0] -= mu_sq * target_left[line];
            rhs_lines[line][ncoef / 2] -= mu_sq * target_right[line];
        }

        for (int i = 0; i < ncoef; ++i) {
            mean_rhs[i] += rhs_lines[line][i];
        }
    }

    for (int i = 0; i < ncoef; ++i) {
        mean_rhs[i] /= Real(nlines);
    }

    auto mean_A = common_A;
    if (!solve_dense_system(mean_A, mean_rhs)) {
        return false;
    }

    auto shifted_A = common_A;
    if (mean_regularization_weight > 0.0) {
        const Real gamma = mean_regularization_weight * scale;
        for (int i = 0; i < ncoef; ++i) {
            shifted_A[i][i] += gamma;
        }

        for (int line = 0; line < nlines; ++line) {
            auto rhs = rhs_lines[line];
            for (int i = 0; i < ncoef; ++i) {
                rhs[i] += gamma * mean_rhs[i];
            }
            auto A = shifted_A;
            if (!solve_dense_system(A, rhs)) {
                return false;
            }
            for (int i = 0; i < ncoef; ++i) {
                alpha_flat[line * ncoef + i] = rhs[i];
            }
        }
    } else {
        for (int line = 0; line < nlines; ++line) {
            auto rhs = rhs_lines[line];
            auto A = common_A;
            if (!solve_dense_system(A, rhs)) {
                return false;
            }
            for (int i = 0; i < ncoef; ++i) {
                alpha_flat[line * ncoef + i] = rhs[i];
            }
        }
    }

    return true;
}

void build_periodic_real_basis(
    int npts,
    int nharmonics,
    std::vector<std::vector<Real>>& basis)
{
    constexpr Real two_pi = 6.28318530717958647692528676655900577;

    const int nmodes = 1 + 2 * nharmonics;
    basis.assign(nmodes, std::vector<Real>(npts, 0.0));

    for (int p = 0; p < npts; ++p) {
        basis[0][p] = 1.0;
        for (int h = 1; h <= nharmonics; ++h) {
            const Real theta = two_pi * Real(h * p) / Real(npts);
            basis[2 * h - 1][p] = std::cos(theta);
            basis[2 * h][p] = std::sin(theta);
        }
    }
}

void project_coefficients_onto_basis(
    std::vector<Real>& coeff_flat,
    int nlines,
    int ncoef,
    const std::vector<std::vector<Real>>& basis)
{
    const int nmodes = static_cast<int>(basis.size());
    if (nmodes == 0 || nlines <= 0 || ncoef <= 0) {
        return;
    }

    std::vector<Real> projected(coeff_flat.size(), 0.0);
    std::vector<Real> gram_diag(nmodes, 0.0);
    for (int m = 0; m < nmodes; ++m) {
        for (int line = 0; line < nlines; ++line) {
            gram_diag[m] += basis[m][line] * basis[m][line];
        }
        gram_diag[m] = std::max(gram_diag[m], Real(1.0e-14));
    }

    for (int r = 0; r < ncoef; ++r) {
        std::vector<Real> modal_coeff(nmodes, 0.0);
        for (int m = 0; m < nmodes; ++m) {
            Real accum = 0.0;
            for (int line = 0; line < nlines; ++line) {
                accum += basis[m][line] * coeff_flat[line * ncoef + r];
            }
            modal_coeff[m] = accum / gram_diag[m];
        }
        for (int line = 0; line < nlines; ++line) {
            Real value = 0.0;
            for (int m = 0; m < nmodes; ++m) {
                value += modal_coeff[m] * basis[m][line];
            }
            projected[line * ncoef + r] = value;
        }
    }

    coeff_flat.swap(projected);
}

bool solve_modal_edge_system(
    const std::vector<std::vector<Real>>& B,
    const std::vector<std::array<Real, 4>>& deltas,
    const std::vector<Real>& target_left,
    const std::vector<Real>& target_right,
    const std::vector<std::vector<Real>>& transverse_basis,
    int ncoef,
    Real boundary_penalty_weight,
    Real coefficient_regularization,
    std::vector<Real>& alpha_flat)
{
    const int nlines = static_cast<int>(deltas.size());
    const int nmodes = static_cast<int>(transverse_basis.size());
    const int nunk = ncoef * nmodes;

    AMREX_ALWAYS_ASSERT(nlines == static_cast<int>(target_left.size()));
    AMREX_ALWAYS_ASSERT(nlines == static_cast<int>(target_right.size()));

    std::vector<std::vector<Real>> M(nunk, std::vector<Real>(nunk, 0.0));
    std::vector<Real> rhs(nunk, 0.0);

    for (int line = 0; line < nlines; ++line) {
        for (int deriv = 0; deriv < 4; ++deriv) {
            for (int r1 = 0; r1 < ncoef; ++r1) {
                for (int m1 = 0; m1 < nmodes; ++m1) {
                    const int idx1 = r1 * nmodes + m1;
                    const Real a1 = transverse_basis[m1][line] * B[deriv][r1];
                    rhs[idx1] -= a1 * deltas[line][deriv];
                    for (int r2 = 0; r2 < ncoef; ++r2) {
                        for (int m2 = 0; m2 < nmodes; ++m2) {
                            const int idx2 = r2 * nmodes + m2;
                            const Real a2 = transverse_basis[m2][line] * B[deriv][r2];
                            M[idx1][idx2] += a1 * a2;
                        }
                    }
                }
            }
        }
    }

    Real scale = 0.0;
    for (int i = 0; i < nunk; ++i) {
        scale += M[i][i];
    }
    scale = (scale / Real(std::max(nunk, 1))) + 1.0e-14;

    const Real lambda_sq = (1.0e-8 + coefficient_regularization) * scale;
    for (int i = 0; i < nunk; ++i) {
        M[i][i] += lambda_sq;
    }

    if (boundary_penalty_weight > 0.0) {
        const Real mu_sq = boundary_penalty_weight * scale;
        const int left_r = 0;
        const int right_r = ncoef / 2;
        for (int line = 0; line < nlines; ++line) {
            for (int m1 = 0; m1 < nmodes; ++m1) {
                const Real b1 = transverse_basis[m1][line];
                const int left_idx1 = left_r * nmodes + m1;
                const int right_idx1 = right_r * nmodes + m1;
                rhs[left_idx1] += mu_sq * b1 * target_left[line];
                rhs[right_idx1] += mu_sq * b1 * target_right[line];
                for (int m2 = 0; m2 < nmodes; ++m2) {
                    const Real b2 = transverse_basis[m2][line];
                    M[left_idx1][left_r * nmodes + m2] += mu_sq * b1 * b2;
                    M[right_idx1][right_r * nmodes + m2] += mu_sq * b1 * b2;
                }
            }
        }
    }

    if (!solve_dense_system(M, rhs)) {
        return false;
    }

    alpha_flat.assign(nlines * ncoef, 0.0);
    for (int line = 0; line < nlines; ++line) {
        for (int r = 0; r < ncoef; ++r) {
            Real alpha = 0.0;
            for (int m = 0; m < nmodes; ++m) {
                alpha += rhs[r * nmodes + m] * transverse_basis[m][line];
            }
            alpha_flat[line * ncoef + r] = alpha;
        }
    }

    return true;
}

} // namespace crescent::nested::c3

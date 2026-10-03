#include "CRESCENT_NestedC3Periodicity.H"

// Active interface-based nested C3 repair implementation.
#include "CRESCENT_NestedC3LinearSolvers.H"

#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <vector>

namespace crescent::nested::c3 {

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real qBasis(Real xi, int r)
{
    const Real one_minus_xi = 1.0 - xi;
    const Real c3_factor =
        one_minus_xi * one_minus_xi * one_minus_xi * one_minus_xi;

    Real xi_power = 1.0;
    for (int m = 0; m < r; ++m) {
        xi_power *= xi;
    }

    return c3_factor * xi_power;
}

struct InterfaceWeights
{
    Real wr[4][5];
    Real wl[4][5];

    InterfaceWeights()
    {
        wr[0][0] =  315.0/128.0;  wr[0][1] = -105.0/32.0;   wr[0][2] =  189.0/64.0;
        wr[0][3] =  -45.0/32.0;   wr[0][4] =   35.0/128.0;

        wr[1][0] =  -31.0/8.0;    wr[1][1] =  229.0/24.0;   wr[1][2] =  -75.0/8.0;
        wr[1][3] =   37.0/8.0;    wr[1][4] =  -11.0/12.0;

        wr[2][0] =  103.0/24.0;   wr[2][1] =  -41.0/3.0;    wr[2][2] =   65.0/4.0;
        wr[2][3] =  -26.0/3.0;    wr[2][4] =   43.0/24.0;

        wr[3][0] =  -3.0;         wr[3][1] =  11.0;         wr[3][2] = -15.0;
        wr[3][3] =   9.0;         wr[3][4] =  -2.0;

        for (int ord = 0; ord < 4; ++ord) {
            const Real sign = (ord % 2 == 0) ? 1.0 : -1.0;
            for (int m = 0; m < 5; ++m) {
                wl[ord][m] = sign * wr[ord][m];
            }
        }
    }
};

const InterfaceWeights g_interface_weights;

void validateInterfaceWeights()
{
    Real sum_r = 0.0;
    Real sum_l = 0.0;
    for (int m = 0; m < 5; ++m) {
        sum_r += g_interface_weights.wr[0][m];
        sum_l += g_interface_weights.wl[0][m];
    }
    if (std::abs(sum_r - 1.0) > 1.0e-10 ||
        std::abs(sum_l - 1.0) > 1.0e-10) {
        amrex::Print() << "ERROR: Interface weight sums incorrect! sum_r="
                       << sum_r << ", sum_l=" << sum_l << std::endl;
        amrex::Abort("Interface weights failed validation");
    }
}

void compute_interface_derivatives_1d(
    const Real* f,
    int n,
    Real dx,
    Real* derivL,
    Real* derivR)
{
    const Real dx_pow[4] = {1.0, dx, dx * dx, dx * dx * dx};

    for (int ord = 0; ord < 4; ++ord) {
        derivR[ord] = 0.0;
        for (int m = 0; m < 5; ++m) {
            derivR[ord] += g_interface_weights.wr[ord][m] * f[m];
        }
        derivR[ord] /= dx_pow[ord];

        derivL[ord] = 0.0;
        for (int m = 0; m < 5; ++m) {
            derivL[ord] += g_interface_weights.wl[ord][m] * f[n - 1 - m];
        }
        derivL[ord] /= dx_pow[ord];
    }
}

const std::vector<std::vector<Real>>& periodicRealBasis(
    crescent::nested::PeriodicRealBasisCache& cache,
    int npts,
    int nharmonics)
{
    if (cache.npts != npts || cache.nharmonics != nharmonics) {
        build_periodic_real_basis(npts, nharmonics, cache.basis);
        cache.npts = npts;
        cache.nharmonics = nharmonics;
    }
    return cache.basis;
}

} // namespace

void computeC3InterfaceBasis(
    crescent::nested::C3RepairWorkspace& workspace,
    int k,
    int nx,
    int ny,
    int bx,
    int by,
    Real dx,
    Real dy)
{
    validateInterfaceWeights();

    workspace.qx_basis.assign(2 * k, std::vector<Real>(nx, 0.0));
    for (int r = 0; r < k; ++r) {
        for (int i = 0; i < nx; ++i) {
            workspace.qx_basis[r][i] =
                (i < bx) ? qBasis(Real(i) / Real(bx), r) : 0.0;
            workspace.qx_basis[k + r][i] =
                (i > nx - 1 - bx)
                    ? qBasis(Real(nx - 1 - i) / Real(bx), r)
                    : 0.0;
        }
    }

    workspace.qy_basis.assign(2 * k, std::vector<Real>(ny, 0.0));
    for (int r = 0; r < k; ++r) {
        for (int j = 0; j < ny; ++j) {
            workspace.qy_basis[r][j] =
                (j < by) ? qBasis(Real(j) / Real(by), r) : 0.0;
            workspace.qy_basis[k + r][j] =
                (j > ny - 1 - by)
                    ? qBasis(Real(ny - 1 - j) / Real(by), r)
                    : 0.0;
        }
    }

    workspace.bx_basis.assign(4, std::vector<Real>(2 * k, 0.0));
    for (int r = 0; r < 2 * k; ++r) {
        Real DL[4], DR[4];
        compute_interface_derivatives_1d(
            workspace.qx_basis[r].data(), nx, dx, DL, DR);
        for (int n = 0; n < 4; ++n) {
            workspace.bx_basis[n][r] = DR[n] - DL[n];
        }
    }

    workspace.by_basis.assign(4, std::vector<Real>(2 * k, 0.0));
    for (int r = 0; r < 2 * k; ++r) {
        Real DL[4], DR[4];
        compute_interface_derivatives_1d(
            workspace.qy_basis[r].data(), ny, dy, DL, DR);
        for (int n = 0; n < 4; ++n) {
            workspace.by_basis[n][r] = DR[n] - DL[n];
        }
    }

    workspace.basis_k = k;
    workspace.basis_nx = nx;
    workspace.basis_ny = ny;
    workspace.basis_bx = bx;
    workspace.basis_by = by;
    workspace.basis_dx = dx;
    workspace.basis_dy = dy;
}

void applyC3InterfaceSeamRepair(
    MultiFab& field,
    const Geometry& geom,
    const BoxArray& ba,
    const C3InterfaceRepairOptions& options,
    crescent::nested::C3RepairWorkspace& workspace)
{
    applyC3InterfaceSeamRepair(
        field,
        geom,
        ba,
        options.buffer_width,
        options.verbosity,
        options.apply_x,
        options.apply_y,
        options.boundary_penalty_weight,
        options.coefficient_regularization_weight,
        options.x_pass_y_smoothing_passes,
        options.x_pass_corner_taper,
        options.x_pass_alpha_y_smoothing_passes,
        options.x_alpha_y_regularization_weight,
        options.y_alpha_x_regularization_weight,
        options.x_alpha_y_mean_regularization_weight,
        options.y_alpha_x_mean_regularization_weight,
        options.x_transverse_harmonics,
        options.y_transverse_harmonics,
        options.x_alpha_y_projection_harmonics,
        options.y_alpha_x_projection_harmonics,
        workspace);
}

void applyC3InterfaceSeamRepair(
    MultiFab& field,
    const Geometry& geom,
    const BoxArray& ba,
    int buffer_width,
    int verbosity,
    bool apply_x,
    bool apply_y,
    Real boundary_penalty_weight,
    Real coefficient_regularization_weight,
    int x_pass_y_smoothing_passes,
    int x_pass_corner_taper,
    int x_pass_alpha_y_smoothing_passes,
    Real x_alpha_y_regularization_weight,
    Real y_alpha_x_regularization_weight,
    Real x_alpha_y_mean_regularization_weight,
    Real y_alpha_x_mean_regularization_weight,
    int x_transverse_harmonics,
    int y_transverse_harmonics,
    int x_alpha_y_projection_harmonics,
    int y_alpha_x_projection_harmonics,
    crescent::nested::C3RepairWorkspace& workspace)
{
    validateInterfaceWeights();

    const int k = InterfaceBasisOrder;
    const Box global_box = ba.minimalBox();
    const int nx = global_box.length(0);
    const int ny = global_box.length(1);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);

    AMREX_ALWAYS_ASSERT(nx >= 5 && ny >= 5);

    const int bx = buffer_width;
    const int by = buffer_width;
    AMREX_ALWAYS_ASSERT(bx >= 8 && by >= 8);

    if (verbosity > 0) {
        amrex::Print() << "\n=== C³-Repair (Interface-based, multi-box) ===="
                       << std::endl;
        amrex::Print() << "  Global box: " << global_box << std::endl;
        amrex::Print() << "  Grid: " << nx << " x " << ny << std::endl;
        amrex::Print() << "  Number of boxes: " << ba.size() << std::endl;
        amrex::Print() << "  dx=" << dx << ", dy=" << dy << std::endl;
        amrex::Print() << "  Buffer width: bx=" << bx << ", by=" << by
                       << std::endl;
        amrex::Print() << "  Basis functions: k=" << k << std::endl;
        amrex::Print() << "  Using 5-point interface reconstruction"
                       << std::endl;
    }

    // Basis matrices Qx, Qy, Bx, By are deterministic in (k, nx, ny, bx, by,
    // dx, dy) and precomputed once at setup. Assert the cached basis still
    // matches this call's parameters; if a future regrid changes any of these
    // the assert trips at the first repair call.
    AMREX_ALWAYS_ASSERT(workspace.basisMatches(k, nx, ny, bx, by, dx, dy));
    const std::vector<std::vector<Real>>& Qx = workspace.qx_basis;
    const std::vector<std::vector<Real>>& Qy = workspace.qy_basis;
    const std::vector<std::vector<Real>>& Bx = workspace.bx_basis;
    const std::vector<std::vector<Real>>& By = workspace.by_basis;
    (void)Qy;  // y-pass reads via workspace.qy_basis further below

    if (verbosity > 0) {
        amrex::Print() << "\n  B_x matrix (2k=" << 2 * k
                       << " basis functions):" << std::endl;
        for (int n = 0; n < 4; ++n) {
            amrex::Print() << "    B_x[" << n << "]: ";
            for (int r = 0; r < 2 * k; ++r) {
                amrex::Print() << Bx[n][r];
                if (r < 2 * k - 1) {
                    amrex::Print() << ", ";
                }
            }
            amrex::Print() << std::endl;
        }
    }

    // Reuse pre-allocated single-box gather target; fixed layout set at setup.
    AMREX_ALWAYS_ASSERT(workspace.ba_global.minimalBox() == global_box);
    MultiFab& field_global = workspace.field_global;
    field_global.ParallelCopy(field, 0, 0, 1);

    for (MFIter mfi(field_global); mfi.isValid(); ++mfi) {
        const Box& bx_box = mfi.validbox();
        auto const& field_arr = field_global.array(mfi);

        AMREX_ALWAYS_ASSERT(bx_box == global_box);

        std::vector<Real>& f_2d = workspace.f_2d_scratch;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int idx = j * nx + i;
                f_2d[idx] = field_arr(
                    i + bx_box.smallEnd(0),
                    j + bx_box.smallEnd(1),
                    0);
            }
        }

        Real max_interface_jump[4] = {0, 0, 0, 0};
        Real max_inner_jump = 0.0;
        const Real coef_reg = coefficient_regularization_weight;

        auto smoothstep01 = [](Real t) -> Real {
            t = std::max(Real(0.0), std::min(Real(1.0), t));
            return t * t * (3.0 - 2.0 * t);
        };

        if (apply_x) {
            std::vector<Real>& x_correction = workspace.x_correction;
            std::vector<Real>& x_alpha = workspace.x_alpha;
            std::vector<std::array<Real, 4>>& x_deltas = workspace.x_deltas;
            std::vector<Real>& x_target_left = workspace.x_target_left;
            std::vector<Real>& x_target_right = workspace.x_target_right;
            std::vector<Real>& x_ref = workspace.x_ref;
            std::vector<Real>& x_u_left = workspace.x_u_left;
            std::vector<Real>& x_u_right = workspace.x_u_right;
            std::fill(x_correction.begin(), x_correction.end(), Real(0.0));
            std::fill(x_alpha.begin(), x_alpha.end(), Real(0.0));
            std::vector<Real>& u_line_buffer = workspace.u_line_scratch;
            Real* u_line = u_line_buffer.data();

            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    u_line[i] = f_2d[j * nx + i];
                }

                Real DL[4], DR[4], Delta[4];
                compute_interface_derivatives_1d(u_line, nx, dx, DL, DR);

                for (int n = 0; n < 4; ++n) {
                    Delta[n] = DR[n] - DL[n];
                    x_deltas[j][n] = Delta[n];
                    max_interface_jump[n] =
                        std::max(max_interface_jump[n], std::abs(Delta[n]));
                }

                const Real ref_x = 0.5 * (u_line[bx] + u_line[nx - 1 - bx]);
                const Real u_left = u_line[0];
                const Real u_right = u_line[nx - 1];
                x_ref[j] = ref_x;
                x_u_left[j] = u_left;
                x_u_right[j] = u_right;
                x_target_left[j] = ref_x - u_left;
                x_target_right[j] = ref_x - u_right;
            }

            bool x_success = true;
            if (x_transverse_harmonics > 0 && ny > 1) {
                const auto& y_basis = periodicRealBasis(
                    workspace.x_pass_transverse_basis,
                    ny,
                    x_transverse_harmonics);
                x_success = solve_modal_edge_system(
                    Bx, x_deltas, x_target_left, x_target_right, y_basis,
                    2 * k, boundary_penalty_weight, coef_reg, x_alpha);
            } else if (x_alpha_y_mean_regularization_weight > 0.0 && ny > 1) {
                x_success = solve_mean_coupled_line_system(
                    Bx, x_deltas, x_target_left, x_target_right, 2 * k,
                    boundary_penalty_weight, coef_reg,
                    x_alpha_y_mean_regularization_weight, x_alpha);
            } else if (x_alpha_y_regularization_weight > 0.0 && ny > 1) {
                x_success = solve_coupled_line_system(
                    Bx, x_deltas, x_target_left, x_target_right, 2 * k,
                    boundary_penalty_weight, coef_reg,
                    x_alpha_y_regularization_weight, x_alpha);
            } else {
                for (int j = 0; j < ny; ++j) {
                    Real alpha[8];
                    const Real line_delta[4] = {
                        x_deltas[j][0], x_deltas[j][1],
                        x_deltas[j][2], x_deltas[j][3]
                    };
                    x_success = solve_least_squares(
                        Bx, line_delta, 2 * k, alpha,
                        boundary_penalty_weight,
                        x_target_left[j], x_target_right[j], coef_reg);
                    if (!x_success) {
                        break;
                    }
                    for (int r = 0; r < 2 * k; ++r) {
                        x_alpha[j * (2 * k) + r] = alpha[r];
                    }
                }
            }

            if (!x_success) {
                amrex::Abort("C3 x-pass coefficient solve failed");
            }

            if (x_transverse_harmonics == 0 &&
                x_alpha_y_projection_harmonics > 0 &&
                ny > 1) {
                const auto& y_basis = periodicRealBasis(
                    workspace.x_pass_projection_basis,
                    ny,
                    x_alpha_y_projection_harmonics);
                project_coefficients_onto_basis(x_alpha, ny, 2 * k, y_basis);
            }

            if (verbosity > 0) {
                amrex::Print() << "\n  First line (j=0) X-PASS diagnostics:"
                               << std::endl;
                amrex::Print() << "    Reference value: ref_x = " << x_ref[0]
                               << std::endl;
                amrex::Print() << "    Before: u[0] = " << x_u_left[0]
                               << ", u[nx-1] = " << x_u_right[0]
                               << std::endl;
                amrex::Print() << "    Targets: delta_L = " << x_target_left[0]
                               << ", delta_R = " << x_target_right[0]
                               << std::endl;
                amrex::Print() << "    Delta: [" << x_deltas[0][0] << ", "
                               << x_deltas[0][1] << ", " << x_deltas[0][2]
                               << ", " << x_deltas[0][3] << "]" << std::endl;
                amrex::Print() << "    Alpha (left buffer r=0.." << k - 1
                               << "): [";
                for (int r = 0; r < k; ++r) {
                    amrex::Print() << x_alpha[r];
                    if (r < k - 1) {
                        amrex::Print() << ", ";
                    }
                }
                amrex::Print() << "]" << std::endl;
                amrex::Print() << "    Alpha (right buffer r=" << k << ".."
                               << 2 * k - 1 << "): [";
                for (int r = k; r < 2 * k; ++r) {
                    amrex::Print() << x_alpha[r];
                    if (r < 2 * k - 1) {
                        amrex::Print() << ", ";
                    }
                }
                amrex::Print() << "]" << std::endl;
            }

            if (x_pass_alpha_y_smoothing_passes > 0 && ny > 1) {
                std::vector<Real>& alpha_work = workspace.x_alpha_work;
                std::fill(alpha_work.begin(), alpha_work.end(), Real(0.0));
                for (int pass = 0; pass < x_pass_alpha_y_smoothing_passes; ++pass) {
                    for (int r = 0; r < 2 * k; ++r) {
                        alpha_work[r] =
                            0.75 * x_alpha[r] + 0.25 * x_alpha[(2 * k) + r];
                        for (int j = 1; j < ny - 1; ++j) {
                            alpha_work[j * (2 * k) + r] =
                                0.25 * x_alpha[(j - 1) * (2 * k) + r] +
                                0.50 * x_alpha[j * (2 * k) + r] +
                                0.25 * x_alpha[(j + 1) * (2 * k) + r];
                        }
                        alpha_work[(ny - 1) * (2 * k) + r] =
                            0.25 * x_alpha[(ny - 2) * (2 * k) + r] +
                            0.75 * x_alpha[(ny - 1) * (2 * k) + r];
                    }
                    x_alpha.swap(alpha_work);
                }
            }

            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    Real delta = 0.0;
                    for (int r = 0; r < 2 * k; ++r) {
                        delta += x_alpha[j * (2 * k) + r] * Qx[r][i];
                    }
                    x_correction[j * nx + i] = delta;
                }
            }

            if (x_pass_y_smoothing_passes > 0 && ny > 1) {
                std::vector<Real>& work = workspace.x_correction_work;
                std::fill(work.begin(), work.end(), Real(0.0));
                for (int pass = 0; pass < x_pass_y_smoothing_passes; ++pass) {
                    for (int i = 0; i < nx; ++i) {
                        work[i] = 0.75 * x_correction[i] +
                                  0.25 * x_correction[nx + i];
                        for (int j = 1; j < ny - 1; ++j) {
                            work[j * nx + i] =
                                0.25 * x_correction[(j - 1) * nx + i] +
                                0.50 * x_correction[j * nx + i] +
                                0.25 * x_correction[(j + 1) * nx + i];
                        }
                        work[(ny - 1) * nx + i] =
                            0.25 * x_correction[(ny - 2) * nx + i] +
                            0.75 * x_correction[(ny - 1) * nx + i];
                    }
                    x_correction.swap(work);
                }
            }

            if (x_pass_corner_taper == 1) {
                for (int j = 0; j < ny; ++j) {
                    Real taper = 1.0;
                    if (j < by) {
                        taper = smoothstep01(Real(j) / Real(by));
                    } else if (j > ny - 1 - by) {
                        taper = smoothstep01(Real(ny - 1 - j) / Real(by));
                    }

                    if (taper < 1.0) {
                        for (int i = 0; i < nx; ++i) {
                            x_correction[j * nx + i] *= taper;
                        }
                    }
                }
            }

            Real* u_line_final = u_line_buffer.data();
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    u_line_final[i] = f_2d[j * nx + i] + x_correction[j * nx + i];
                }

                if (j == 0 && verbosity > 0) {
                    amrex::Print() << "    After: u[0] = " << u_line_final[0]
                                   << ", u[nx-1] = " << u_line_final[nx - 1]
                                   << std::endl;
                }

                if (bx < nx / 2) {
                    const Real sL = (u_line_final[bx] - u_line_final[bx - 1]) / dx;
                    const Real sR = (u_line_final[bx + 1] - u_line_final[bx]) / dx;
                    max_inner_jump =
                        std::max(max_inner_jump, std::abs(sR - sL));
                }

                for (int i = 0; i < nx; ++i) {
                    f_2d[j * nx + i] = u_line_final[i];
                }
            }
        }

        if (apply_y) {
            std::vector<Real>& y_alpha = workspace.y_alpha;
            std::vector<std::array<Real, 4>>& y_deltas = workspace.y_deltas;
            std::vector<Real>& y_target_bottom = workspace.y_target_bottom;
            std::vector<Real>& y_target_top = workspace.y_target_top;
            std::fill(y_alpha.begin(), y_alpha.end(), Real(0.0));
            Real* u_line_y = workspace.u_line_scratch.data();

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    u_line_y[j] = f_2d[j * nx + i];
                }

                Real DL[4], DR[4], Delta[4];
                compute_interface_derivatives_1d(u_line_y, ny, dy, DL, DR);
                for (int n = 0; n < 4; ++n) {
                    Delta[n] = DR[n] - DL[n];
                    y_deltas[i][n] = Delta[n];
                }

                const Real ref_y = 0.5 * (u_line_y[by] + u_line_y[ny - 1 - by]);
                y_target_bottom[i] = ref_y - u_line_y[0];
                y_target_top[i] = ref_y - u_line_y[ny - 1];
            }

            bool y_success = true;
            if (y_transverse_harmonics > 0 && nx > 1) {
                const auto& x_basis = periodicRealBasis(
                    workspace.y_pass_transverse_basis,
                    nx,
                    y_transverse_harmonics);
                y_success = solve_modal_edge_system(
                    By, y_deltas, y_target_bottom, y_target_top, x_basis,
                    2 * k, boundary_penalty_weight, coef_reg, y_alpha);
            } else if (y_alpha_x_mean_regularization_weight > 0.0 && nx > 1) {
                y_success = solve_mean_coupled_line_system(
                    By, y_deltas, y_target_bottom, y_target_top, 2 * k,
                    boundary_penalty_weight, coef_reg,
                    y_alpha_x_mean_regularization_weight, y_alpha);
            } else if (y_alpha_x_regularization_weight > 0.0 && nx > 1) {
                y_success = solve_coupled_line_system(
                    By, y_deltas, y_target_bottom, y_target_top, 2 * k,
                    boundary_penalty_weight, coef_reg,
                    y_alpha_x_regularization_weight, y_alpha);
            } else {
                for (int i = 0; i < nx; ++i) {
                    Real alpha[8];
                    const Real line_delta[4] = {
                        y_deltas[i][0], y_deltas[i][1],
                        y_deltas[i][2], y_deltas[i][3]
                    };
                    y_success = solve_least_squares(
                        By, line_delta, 2 * k, alpha,
                        boundary_penalty_weight,
                        y_target_bottom[i], y_target_top[i], coef_reg);
                    if (!y_success) {
                        break;
                    }
                    for (int r = 0; r < 2 * k; ++r) {
                        y_alpha[i * (2 * k) + r] = alpha[r];
                    }
                }
            }

            if (!y_success) {
                amrex::Abort("C3 y-pass coefficient solve failed");
            }

            if (y_transverse_harmonics == 0 &&
                y_alpha_x_projection_harmonics > 0 &&
                nx > 1) {
                const auto& x_basis = periodicRealBasis(
                    workspace.y_pass_projection_basis,
                    nx,
                    y_alpha_x_projection_harmonics);
                project_coefficients_onto_basis(y_alpha, nx, 2 * k, x_basis);
            }

            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    Real delta = 0.0;
                    for (int r = 0; r < 2 * k; ++r) {
                        delta += y_alpha[i * (2 * k) + r] * Qy[r][j];
                    }
                    f_2d[j * nx + i] += delta;
                }
            }
        }

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                field_arr(
                    i + bx_box.smallEnd(0),
                    j + bx_box.smallEnd(1),
                    0) = f_2d[j * nx + i];
            }
        }

        if (verbosity > 0) {
            amrex::Print() << "\n  Interface mismatch (before correction):"
                           << std::endl;
            amrex::Print() << "    |Δ₀| = " << max_interface_jump[0]
                           << std::endl;
            amrex::Print() << "    |Δ₁| = " << max_interface_jump[1]
                           << std::endl;
            amrex::Print() << "    |Δ₂| = " << max_interface_jump[2]
                           << std::endl;
            amrex::Print() << "    |Δ₃| = " << max_interface_jump[3]
                           << std::endl;
            amrex::Print() << "  Inner-interface slope jump: " << max_inner_jump
                           << std::endl;

            Real max_after[4] = {0, 0, 0, 0};
            Real* u_line_after = workspace.u_line_scratch.data();
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    u_line_after[i] = f_2d[j * nx + i];
                }
                Real DL[4], DR[4], Delta[4];
                compute_interface_derivatives_1d(u_line_after, nx, dx, DL, DR);
                for (int n = 0; n < 4; ++n) {
                    Delta[n] = DR[n] - DL[n];
                    max_after[n] = std::max(max_after[n], std::abs(Delta[n]));
                }
            }
            amrex::Print() << "\n  Interface mismatch (after correction):"
                           << std::endl;
            amrex::Print() << "    |Δ₀| = " << max_after[0] << std::endl;
            amrex::Print() << "    |Δ₁| = " << max_after[1] << std::endl;
            amrex::Print() << "    |Δ₂| = " << max_after[2] << std::endl;
            amrex::Print() << "    |Δ₃| = " << max_after[3] << std::endl;
            amrex::Print() << "  Reduction: " << max_interface_jump[0]
                           << " → " << max_after[0] << " ("
                           << (max_interface_jump[0] /
                               std::max(max_after[0], 1.0e-10))
                           << "×)" << std::endl;
        }
    }

    field.ParallelCopy(field_global, 0, 0, 1);

    if (verbosity > 0) {
        amrex::Print()
            << "=== C³-Repair (Interface-based, multi-box) complete ==="
            << std::endl;
    }
}

} // namespace crescent::nested::c3

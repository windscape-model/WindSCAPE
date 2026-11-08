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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

#ifndef AMREX_USE_GPU
// Least-squares coefficients of the 2k basis functions that remove the four
// wrap jumps of one line (value and first three derivatives) and pull both
// line ends toward the mean of the two buffer-edge values.
bool solveLineCoefficients(
    const Real* u,
    int n,
    Real h,
    int b,
    const std::vector<std::vector<Real>>& B,
    int ncoef,
    Real boundary_penalty_weight,
    Real coefficient_regularization,
    Real alpha[])
{
    Real DL[4], DR[4], Delta[4];
    compute_interface_derivatives_1d(u, n, h, DL, DR);
    for (int ord = 0; ord < 4; ++ord) {
        Delta[ord] = DR[ord] - DL[ord];
    }
    const Real ref = 0.5 * (u[b] + u[n - 1 - b]);
    return solve_least_squares(
        B, Delta, ncoef, alpha, boundary_penalty_weight,
        ref - u[0], ref - u[n - 1], coefficient_regularization);
}
#endif

#ifdef AMREX_USE_GPU
// Batched GPU path for the default linewise interface repair. The host path
// below is unchanged and remains the CPU implementation; this runs only on GPU
// builds. Each grid line is independent, so a pass is a ParallelFor over lines.
struct C3Weights { amrex::Real wr[4][5]; amrex::Real wl[4][5]; };
struct C3Basis   { amrex::Real m[4][8]; };  // 4 x 2k basis-jump matrix (2k<=8)

void c3SeamRepairPassDevice(
    amrex::Array4<amrex::Real> const& f, const amrex::Box& gb,
    int dir, int nx, int ny, int b, amrex::Real dl,
    C3Weights w, C3Basis Bmat,
    amrex::Real penalty, amrex::Real coef_reg, int k)
{
    using amrex::Real;
    const int i0 = gb.smallEnd(0), j0 = gb.smallEnd(1), k0 = gb.smallEnd(2);
    const int nlines = (dir == 0) ? ny : nx;   // rows for x-pass, cols for y-pass
    const int n = (dir == 0) ? nx : ny;         // length along the pass direction
    amrex::ParallelFor(nlines, [=] AMREX_GPU_DEVICE(int line) noexcept {
        auto at = [=] (int p) -> Real {
            return (dir == 0) ? f(i0 + p, j0 + line, k0)
                              : f(i0 + line, j0 + p, k0);
        };
        Real Delta[4];
        for (int ord = 0; ord < 4; ++ord) {
            Real rr = 0.0, ll = 0.0;
            for (int m = 0; m < 5; ++m) {
                rr += w.wr[ord][m] * at(m);
                ll += w.wl[ord][m] * at(n - 1 - m);
            }
            Real dp = 1.0;
            for (int q = 0; q < ord; ++q) { dp *= dl; }
            Delta[ord] = (rr - ll) / dp;
        }
        const Real ref = 0.5 * (at(b) + at(n - 1 - b));
        const Real tl = ref - at(0);
        const Real tr = ref - at(n - 1);
        Real alpha[8] = {};
        solve_least_squares_device(Bmat.m, Delta, 2 * k, alpha, penalty, tl, tr,
                                   coef_reg);
        for (int p = 0; p < b; ++p) {
            Real cL = 0.0, cR = 0.0;
            for (int r = 0; r < k; ++r) {
                const Real q = qBasis(Real(p) / Real(b), r);
                cL += alpha[r] * q;
                cR += alpha[k + r] * q;
            }
            if (dir == 0) {
                f(i0 + p,         j0 + line, k0) += cL;
                f(i0 + (n - 1 - p), j0 + line, k0) += cR;
            } else {
                f(i0 + line, j0 + p,           k0) += cL;
                f(i0 + line, j0 + (n - 1 - p), k0) += cR;
            }
        }
    });
}
#endif // AMREX_USE_GPU

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
    const int k = InterfaceBasisOrder;
    const Box global_box = ba.minimalBox();
    const int nx = global_box.length(0);
    const int ny = global_box.length(1);

    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);

    AMREX_ALWAYS_ASSERT(nx >= 5 && ny >= 5);

    const int bx = options.buffer_width_x;
    const int by = options.buffer_width_y;
    AMREX_ALWAYS_ASSERT(bx >= 8 && by >= 8);
    const Real penalty = options.boundary_penalty_weight;
    const Real coef_reg = options.coefficient_regularization_weight;

    // Basis matrices Qx, Qy, Bx, By are deterministic in (k, nx, ny, bx, by,
    // dx, dy) and precomputed once at setup. Assert the cached basis still
    // matches this call's parameters; if a future regrid changes any of these
    // the assert trips at the first repair call.
    AMREX_ALWAYS_ASSERT(workspace.basisMatches(k, nx, ny, bx, by, dx, dy));
    const std::vector<std::vector<Real>>& Qx = workspace.qx_basis;
    const std::vector<std::vector<Real>>& Qy = workspace.qy_basis;
    const std::vector<std::vector<Real>>& Bx = workspace.bx_basis;
    const std::vector<std::vector<Real>>& By = workspace.by_basis;

    // Reuse pre-allocated single-box gather target; fixed layout set at setup.
    AMREX_ALWAYS_ASSERT(workspace.ba_global.minimalBox() == global_box);
    MultiFab& field_global = workspace.field_global;
    field_global.ParallelCopy(field, 0, 0, 1);

#ifdef AMREX_USE_GPU
    {
        C3Weights w;
        for (int o = 0; o < 4; ++o) {
            for (int m = 0; m < 5; ++m) {
                w.wr[o][m] = g_interface_weights.wr[o][m];
                w.wl[o][m] = g_interface_weights.wl[o][m];
            }
        }
        C3Basis Bxd{}, Byd{};
        for (int nn = 0; nn < 4; ++nn) {
            for (int r = 0; r < 2 * k; ++r) {
                Bxd.m[nn][r] = Bx[nn][r];
                Byd.m[nn][r] = By[nn][r];
            }
        }
        for (MFIter mfi(field_global); mfi.isValid(); ++mfi) {
            auto arr = field_global.array(mfi);
            const Box& gb = mfi.validbox();
            c3SeamRepairPassDevice(arr, gb, 0, nx, ny, bx, dx, w, Bxd,
                                   penalty, coef_reg, k);
            c3SeamRepairPassDevice(arr, gb, 1, nx, ny, by, dy, w, Byd,
                                   penalty, coef_reg, k);
        }
    }
#else
    for (MFIter mfi(field_global); mfi.isValid(); ++mfi) {
        const Box& bx_box = mfi.validbox();
        auto const& field_arr = field_global.array(mfi);

        AMREX_ALWAYS_ASSERT(bx_box == global_box);

        std::vector<Real>& f_2d = workspace.f_2d_scratch;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                f_2d[j * nx + i] = field_arr(
                    i + bx_box.smallEnd(0),
                    j + bx_box.smallEnd(1),
                    0);
            }
        }

        Real* u_line = workspace.u_line_scratch.data();
        Real alpha[8];

        // x pass: each row is fitted and corrected independently.
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                u_line[i] = f_2d[j * nx + i];
            }
            if (!solveLineCoefficients(u_line, nx, dx, bx, Bx, 2 * k,
                                       penalty, coef_reg, alpha)) {
                amrex::Abort("C3 x-pass coefficient solve failed");
            }
            for (int i = 0; i < nx; ++i) {
                Real delta = 0.0;
                for (int r = 0; r < 2 * k; ++r) {
                    delta += alpha[r] * Qx[r][i];
                }
                f_2d[j * nx + i] = f_2d[j * nx + i] + delta;
            }
        }

        // y pass on the x-repaired field: each column independently.
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                u_line[j] = f_2d[j * nx + i];
            }
            if (!solveLineCoefficients(u_line, ny, dy, by, By, 2 * k,
                                       penalty, coef_reg, alpha)) {
                amrex::Abort("C3 y-pass coefficient solve failed");
            }
            for (int j = 0; j < ny; ++j) {
                Real delta = 0.0;
                for (int r = 0; r < 2 * k; ++r) {
                    delta += alpha[r] * Qy[r][j];
                }
                f_2d[j * nx + i] += delta;
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
    }
#endif // AMREX_USE_GPU

    field.ParallelCopy(field_global, 0, 0, 1);
}

} // namespace crescent::nested::c3

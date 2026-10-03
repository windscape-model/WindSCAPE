#include "CRESCENT_NestedMasks.H"
#include <cmath>

namespace crescent::nested::masks {

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real septic_smoothstep(Real t)
{
    t = amrex::max(Real(0.0), amrex::min(Real(1.0), t));
    const Real t2 = t * t;
    const Real t4 = t2 * t2;
    return t4 * (35.0 - 84.0 * t + 70.0 * t2 - 20.0 * t2 * t);
}

Real initial_euclidean_sigma(const domains::DomainStructure& domain, Real x, Real y)
{
    Real dist_x = 0.0;
    if (x < domain.x_lo_A) {
        dist_x = domain.x_lo_A - x;
    } else if (x > domain.x_hi_A) {
        dist_x = x - domain.x_hi_A;
    }

    Real dist_y = 0.0;
    if (y < domain.y_lo_A) {
        dist_y = domain.y_lo_A - y;
    } else if (y > domain.y_hi_A) {
        dist_y = y - domain.y_hi_A;
    }

    if (dist_x == 0.0 && dist_y == 0.0) {
        return 0.0;
    }

    Real normalized_dist = 0.0;
    if (dist_x > 0.0 && dist_y > 0.0) {
        const Real nx = dist_x / domain.W_mix_x;
        const Real ny = dist_y / domain.W_mix_y;
        normalized_dist = std::sqrt(nx * nx + ny * ny);
    } else if (dist_x > 0.0) {
        normalized_dist = dist_x / domain.W_mix_x;
    } else {
        normalized_dist = dist_y / domain.W_mix_y;
    }

    return (normalized_dist >= 1.0) ? 1.0 : septic_smoothstep(normalized_dist);
}

Real blending_ramp_1d(Real coord, Real a_lo, Real a_hi, Real s_lo, Real s_hi)
{
    if (coord < a_lo) {
        const Real width = a_lo - s_lo;
        if (width <= 0.0) {
            return 1.0;
        }
        const Real normalized_dist = (a_lo - coord) / width;
        return (normalized_dist >= 1.0) ? 1.0 : septic_smoothstep(normalized_dist);
    }

    if (coord > a_hi) {
        const Real width = s_hi - a_hi;
        if (width <= 0.0) {
            return 1.0;
        }
        const Real normalized_dist = (coord - a_hi) / width;
        return (normalized_dist >= 1.0) ? 1.0 : septic_smoothstep(normalized_dist);
    }

    return 0.0;
}

Real blending_transition_coordinate_1d(Real coord, Real a_lo, Real a_hi, Real s_lo, Real s_hi)
{
    if (coord < a_lo) {
        const Real width = a_lo - s_lo;
        if (width <= 0.0) {
            return 1.0;
        }
        const Real normalized_dist = (a_lo - coord) / width;
        return amrex::min(Real(1.0), amrex::max(Real(0.0), normalized_dist));
    }

    if (coord > a_hi) {
        const Real width = s_hi - a_hi;
        if (width <= 0.0) {
            return 1.0;
        }
        const Real normalized_dist = (coord - a_hi) / width;
        return amrex::min(Real(1.0), amrex::max(Real(0.0), normalized_dist));
    }

    return 0.0;
}

void solve_harmonic_sigma_field(
    std::vector<Real>& sigma,
    std::vector<int>& fixed,
    const domains::DomainStructure& domain,
    const Box& global_box,
    const Geometry& geom)
{
    const int i0 = global_box.smallEnd(0);
    const int j0 = global_box.smallEnd(1);
    const int nx = global_box.length(0);
    const int ny = global_box.length(1);
    const Real* dx = geom.CellSize();
    const auto& prob_lo = geom.ProbLoArray();

    auto idx = [nx](int ii, int jj) { return jj * nx + ii; };

    const int a_lo_i = domain.box_A.smallEnd(0);
    const int a_hi_i = domain.box_A.bigEnd(0);
    const int a_lo_j = domain.box_A.smallEnd(1);
    const int a_hi_j = domain.box_A.bigEnd(1);
    const int s_lo_i = domain.box_S.smallEnd(0);
    const int s_hi_i = domain.box_S.bigEnd(0);
    const int s_lo_j = domain.box_S.smallEnd(1);
    const int s_hi_j = domain.box_S.bigEnd(1);

    int n_unknown = 0;
    for (int j = 0; j < ny; ++j) {
        const int jg = j0 + j;
        for (int i = 0; i < nx; ++i) {
            const int ig = i0 + i;
            const int flat = idx(i, j);

            const bool in_A = (ig >= a_lo_i && ig <= a_hi_i &&
                               jg >= a_lo_j && jg <= a_hi_j);
            const bool on_or_outside_S = (ig <= s_lo_i || ig >= s_hi_i ||
                                          jg <= s_lo_j || jg >= s_hi_j);

            if (in_A) {
                sigma[flat] = 0.0;
                fixed[flat] = 1;
            } else if (on_or_outside_S) {
                sigma[flat] = 1.0;
                fixed[flat] = 1;
            } else {
                const Real x = prob_lo[0] + (static_cast<Real>(ig) + 0.5) * dx[0];
                const Real y = prob_lo[1] + (static_cast<Real>(jg) + 0.5) * dx[1];
                sigma[flat] = initial_euclidean_sigma(domain, x, y);
                fixed[flat] = 0;
                ++n_unknown;
            }
        }
    }

    if (n_unknown <= 0) {
        return;
    }

    const Real inv_dx2 = 1.0 / (dx[0] * dx[0]);
    const Real inv_dy2 = 1.0 / (dx[1] * dx[1]);
    const Real denom = 2.0 * (inv_dx2 + inv_dy2);
    const Real omega = 1.85;
    const int max_iters = 20000;
    const Real tol = 1.0e-10;

    for (int iter = 0; iter < max_iters; ++iter) {
        Real max_delta = 0.0;

        for (int color = 0; color < 2; ++color) {
            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {
                    if (((i + j) & 1) != color) {
                        continue;
                    }

                    const int flat = idx(i, j);
                    if (fixed[flat] != 0) {
                        continue;
                    }

                    const Real harmonic =
                        ((sigma[idx(i - 1, j)] + sigma[idx(i + 1, j)]) * inv_dx2 +
                         (sigma[idx(i, j - 1)] + sigma[idx(i, j + 1)]) * inv_dy2) / denom;
                    const Real updated = amrex::max(Real(0.0),
                                                    amrex::min(Real(1.0),
                                                               sigma[flat] + omega * (harmonic - sigma[flat])));
                    max_delta = std::max(max_delta, std::abs(updated - sigma[flat]));
                    sigma[flat] = updated;
                }
            }
        }

        if (max_delta < tol) {
            break;
        }
    }
}

} // namespace

void generateBlendingFunction(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain,
    Real steepness)
{
    (void) steepness;

    const auto& prob_lo = geom.ProbLoArray();
    const Real* dx = geom.CellSize();

    for (MFIter mfi(sigma_blend); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto sigma = sigma_blend.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];

            // Distance-based blending approach to avoid corner artifacts
            // Compute distance components from point to Ω_A boundary

            Real dist_x = 0.0;
            if (x < domain.x_lo_A) {
                dist_x = domain.x_lo_A - x;  // Distance to left edge of Ω_A
            } else if (x > domain.x_hi_A) {
                dist_x = x - domain.x_hi_A;  // Distance to right edge of Ω_A
            }
            // else: x is inside Ω_A in x-direction, dist_x = 0

            Real dist_y = 0.0;
            if (y < domain.y_lo_A) {
                dist_y = domain.y_lo_A - y;  // Distance to bottom edge of Ω_A
            } else if (y > domain.y_hi_A) {
                dist_y = y - domain.y_hi_A;  // Distance to top edge of Ω_A
            }
            // else: y is inside Ω_A in y-direction, dist_y = 0

            // Inside Ω_A: sigma = 0 (pure L1 spectral solution)
            if (dist_x == 0.0 && dist_y == 0.0) {
                sigma(i,j,k) = 0.0;
            } else {
                // Outside Ω_A: compute normalized distance for blending
                // Normalize by mixing zone widths to account for anisotropy

                Real normalized_dist;
                if (dist_x > 0.0 && dist_y > 0.0) {
                    // Corner region: use Euclidean distance in normalized space
                    // This gives isotropic blending in corner regions
                    Real nx = dist_x / domain.W_mix_x;
                    Real ny = dist_y / domain.W_mix_y;
                    normalized_dist = std::sqrt(nx * nx + ny * ny);
                } else if (dist_x > 0.0) {
                    // Left/right side region: only x-distance matters
                    normalized_dist = dist_x / domain.W_mix_x;
                } else {
                    // Top/bottom side region: only y-distance matters
                    normalized_dist = dist_y / domain.W_mix_y;
                }

                // Check if outside Ω_S (exterior region): pure L0
                bool outside_S = (x < domain.x_lo_S) || (x > domain.x_hi_S) ||
                                 (y < domain.y_lo_S) || (y > domain.y_hi_S);

                if (outside_S) {
                    sigma(i,j,k) = 1.0;
                } else if (normalized_dist >= 1.0) {
                    // At or beyond Ω_S edge within the mixing zone
                    sigma(i,j,k) = 1.0;
                } else {
                    // Inside mixing zone: smooth transition from 0 at Ω_A edge to 1 at Ω_S edge.
                    sigma(i,j,k) = septic_smoothstep(normalized_dist);
                }
            }
        });
    }
}

void generateBlendingFunctionSequential(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain,
    Real steepness)
{
    (void) steepness;

    const auto& prob_lo = geom.ProbLoArray();
    const Real* dx = geom.CellSize();

    for (MFIter mfi(sigma_blend); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto sigma = sigma_blend.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];

            // Sequential 1D blending: compute σ_x and σ_y independently,
            // then take the maximum to ensure proper blending at boundaries

            // ===================================================================
            // X-direction blending (for all y, including corners)
            // ===================================================================
            Real sigma_x = 0.0;
            Real dist_x = 0.0;

            if (x < domain.x_lo_A) {
                dist_x = domain.x_lo_A - x;
            } else if (x > domain.x_hi_A) {
                dist_x = x - domain.x_hi_A;
            }

            if (dist_x > 0.0) {
                Real normalized_dist_x = dist_x / domain.W_mix_x;

                // Check if outside Ω_S in x-direction
                if (x < domain.x_lo_S || x > domain.x_hi_S) {
                    sigma_x = 1.0;
                } else if (normalized_dist_x >= 1.0) {
                    sigma_x = 1.0;
                } else {
                    sigma_x = septic_smoothstep(normalized_dist_x);
                }
            }

            // ===================================================================
            // Y-direction blending (for all x, including corners)
            // ===================================================================
            Real sigma_y = 0.0;
            Real dist_y = 0.0;

            if (y < domain.y_lo_A) {
                dist_y = domain.y_lo_A - y;
            } else if (y > domain.y_hi_A) {
                dist_y = y - domain.y_hi_A;
            }

            if (dist_y > 0.0) {
                Real normalized_dist_y = dist_y / domain.W_mix_y;

                // Check if outside Ω_S in y-direction
                if (y < domain.y_lo_S || y > domain.y_hi_S) {
                    sigma_y = 1.0;
                } else if (normalized_dist_y >= 1.0) {
                    sigma_y = 1.0;
                } else {
                    sigma_y = septic_smoothstep(normalized_dist_y);
                }
            }

            // ===================================================================
            // Combine with smooth transition at boundaries
            // ===================================================================
            // Strategy: Use priority logic with smooth transition zones to reduce
            // discontinuities at Ω_A boundaries where control switches between
            // x-buffer and y-buffer.

            // Define transition zone width (in units of grid cells)
            // Use ~2-3 cells to smooth the transition
            Real transition_width_x = 3.0 * dx[0];
            Real transition_width_y = 3.0 * dx[1];

            if (dist_x > 0.0) {
                // ===== X-BUFFER ZONE =====
                // X is in buffer → use σ_x (for ALL y, including corners) → y-uniform
                sigma(i,j,k) = sigma_x;

            } else if (dist_x == 0.0 && dist_y > 0.0) {
                // ===== TRANSITION AT Ω_A X-BOUNDARY =====
                // X just crossed into Ω_A, but Y is in buffer
                // Smoothly transition from σ_x control to σ_y control

                // Distance from x to nearest x-buffer
                Real dist_to_x_buffer = std::min(x - domain.x_lo_A, domain.x_hi_A - x);

                if (dist_to_x_buffer < transition_width_x) {
                    // In transition zone: blend σ_x (extended) with σ_y
                    // Compute σ_x as if we're just inside the buffer
                    Real virtual_dist_x = transition_width_x - dist_to_x_buffer;
                    Real normalized_virt_x = virtual_dist_x / domain.W_mix_x;
                    Real virtual_sigma_x = 0.0;

                    if (normalized_virt_x < 1.0) {
                        virtual_sigma_x = septic_smoothstep(normalized_virt_x);
                    }

                    // Weight: favors σ_x near boundary, σ_y deeper inside Ω_A
                    Real weight_x = dist_to_x_buffer / transition_width_x;  // 0 at boundary, 1 deeper inside
                    weight_x = weight_x * weight_x;  // Square for smoother transition

                    // Blend: near boundary use extended σ_x, deeper inside use σ_y
                    sigma(i,j,k) = (1.0 - weight_x) * virtual_sigma_x + weight_x * sigma_y;

                } else {
                    // Far from x-buffer: pure σ_y control
                    sigma(i,j,k) = sigma_y;
                }

            } else {
                // ===== INSIDE Ω_A (both x and y) =====
                // Pure L1 spectral solution
                sigma(i,j,k) = 0.0;
            }
        });
    }
}

void generateBlendingFunctionHarmonic(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain,
    Real steepness)
{
    (void) steepness;

    const Box& global_box = sigma_blend.boxArray().minimalBox();
    const int i0 = global_box.smallEnd(0);
    const int j0 = global_box.smallEnd(1);
    const int nx = global_box.length(0);
    const int ny = global_box.length(1);

    std::vector<Real> sigma(nx * ny, 1.0);
    std::vector<int> fixed(nx * ny, 1);
    solve_harmonic_sigma_field(sigma, fixed, domain, global_box, geom);

    for (MFIter mfi(sigma_blend); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto sigma_arr = sigma_blend.array(mfi);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jj = j - j0;
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ii = i - i0;
                sigma_arr(i, j, 0) = sigma[jj * nx + ii];
            }
        }
    }
}

void generateBlendingFunctionHarmonicSmooth(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain,
    Real steepness)
{
    (void) steepness;

    const Box& global_box = sigma_blend.boxArray().minimalBox();
    const int i0 = global_box.smallEnd(0);
    const int j0 = global_box.smallEnd(1);
    const int nx = global_box.length(0);
    const int ny = global_box.length(1);

    std::vector<Real> sigma(nx * ny, 1.0);
    std::vector<int> fixed(nx * ny, 1);
    solve_harmonic_sigma_field(sigma, fixed, domain, global_box, geom);

    for (Real& value : sigma) {
        value = septic_smoothstep(value);
    }

    for (MFIter mfi(sigma_blend); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto sigma_arr = sigma_blend.array(mfi);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jj = j - j0;
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ii = i - i0;
                sigma_arr(i, j, 0) = sigma[jj * nx + ii];
            }
        }
    }
}

void generateBlendingFunctionTensorProduct(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain,
    Real alpha)
{
    const Real safe_alpha = amrex::max(alpha, Real(1.0e-12));

    const auto& prob_lo = geom.ProbLoArray();
    const Real* dx = geom.CellSize();

    for (MFIter mfi(sigma_blend); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto sigma = sigma_blend.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];

            const Real sigma_x = blending_ramp_1d(x, domain.x_lo_A, domain.x_hi_A,
                                                  domain.x_lo_S, domain.x_hi_S);
            const Real sigma_y = blending_ramp_1d(y, domain.y_lo_A, domain.y_hi_A,
                                                  domain.y_lo_S, domain.y_hi_S);

            // Inclusion-exclusion composition preserves the endpoint constraints
            // σ=0 on Ω_A and σ=1 on Ω_S/Ω_E. For alpha=1 this reduces to the
            // standard tensor-product form. For alpha>1 the profile stays closer
            // to zero near Ω_A and turns on later in both the side ramps and
            // the corner regions.
            const Real one_minus_sigma_x = amrex::max(Real(0.0), Real(1.0) - sigma_x);
            const Real one_minus_sigma_y = amrex::max(Real(0.0), Real(1.0) - sigma_y);
            sigma(i, j, k) = 1.0
                - std::pow(one_minus_sigma_x, safe_alpha)
                  * std::pow(one_minus_sigma_y, safe_alpha);
        });
    }
}

void generateBlendingFunctionSuperellipse(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain,
    Real gamma,
    Real p)
{
    const Real safe_gamma = amrex::max(gamma, Real(1.0e-12));
    const Real safe_p = amrex::max(p, Real(1.0e-12));

    const auto& prob_lo = geom.ProbLoArray();
    const Real* dx = geom.CellSize();

    for (MFIter mfi(sigma_blend); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto sigma = sigma_blend.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const Real y = prob_lo[1] + (j + 0.5) * dx[1];

            const Real tx = blending_transition_coordinate_1d(
                x, domain.x_lo_A, domain.x_hi_A, domain.x_lo_S, domain.x_hi_S);
            const Real ty = blending_transition_coordinate_1d(
                y, domain.y_lo_A, domain.y_hi_A, domain.y_lo_S, domain.y_hi_S);

            const Real ux = std::pow(tx, safe_gamma);
            const Real uy = std::pow(ty, safe_gamma);
            const Real radius =
                std::pow(std::pow(ux, safe_p) + std::pow(uy, safe_p), 1.0 / safe_p);
            sigma(i, j, k) = septic_smoothstep(amrex::min(Real(1.0), radius));
        });
    }
}

} // namespace crescent::nested::masks

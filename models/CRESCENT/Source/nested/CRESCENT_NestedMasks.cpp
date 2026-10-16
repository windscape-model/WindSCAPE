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

// One-dimensional ramp: 0 inside [a_lo, a_hi], the septic smoothstep of the
// normalized distance across the mixing annulus, and 1 at or beyond the
// Omega_S boundary s_lo / s_hi.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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

} // namespace

void generateBlendingFunction(
    MultiFab& sigma_blend,
    const Geometry& geom,
    const domains::DomainStructure& domain)
{
    const auto& prob_lo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

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

            const Real one_minus_sigma_x = amrex::max(Real(0.0), Real(1.0) - sigma_x);
            const Real one_minus_sigma_y = amrex::max(Real(0.0), Real(1.0) - sigma_y);
            sigma(i, j, k) = 1.0 - one_minus_sigma_x * one_minus_sigma_y;
        });
    }
}

} // namespace crescent::nested::masks

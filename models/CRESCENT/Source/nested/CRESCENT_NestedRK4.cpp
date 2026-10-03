#include "CRESCENT_NestedRK4.H"

// Active nested-L1 RK4 implementation.
#include "CRESCENT_Interpolation.H"
#include "CRESCENT_Spectral.H"
#include <AMReX_Print.H>
#include <cmath>
#include <utility>

namespace crescent::nested::rk4 {

using namespace amrex;

// ============================================================================
// Spatial + temporal interpolation of L0 derivatives to L1 for RK4
// ============================================================================
// Purpose: Interpolate L0 derivatives to L1 grid at RK4 substage time
// Usage: During CRESCENT timestepping (called 5× per RK4 stage for each derivative)
//
// Key features:
// - Spatial interpolation: L0 grid → L1 grid (selectable: bilinear or quadratic)
// - Temporal interpolation: Blend between t^n and t^{n+1}
//   Formula: D(t^n + c·Δt) = (1-c)·D^n + c·D^{n+1}
//   where c = RK4 coefficient [0.0, 0.5, 0.5, 1.0] for stages 0-3
// - Reuses the custom prolongation path used elsewhere in the codebase
// - nz=1 compatible (2D physics in 3D build)
//
// Multi-box safety:
// - Before fix: Directly accessed L0 boxes → some L1 regions got wrong data
// - After fix: ParallelCopy gathers all L0 data into single patch first
// - Strategy matches crescent::interpolation::interpolateParentToChild()
//
// Related function:
// - crescent::interpolation::interpolateParentToChild() handles spatial-only
//   interpolation for state variables (eta, phi) at initialization.
//
// Why separate functions:
// - This function: spatial + temporal (needs two time levels for RK4)
// - Initialization helper: Spatial only (single time level, initialization)
// ============================================================================
void interpolateCoreDerivativeToL1Stage(
    const MultiFab& deriv_L0_n,
    const MultiFab& deriv_L0_np1,
    Real c,
    const Geometry& geom_L0,
    const Geometry& geom_L1,
    MultiFab& deriv_interp,
    MultiFab& deriv_stage_scratch,
    MultiFab& coarse_patch_scratch,
    int interpolation_order)
{
    // Temporal interpolation coefficient
    // t_substage = t^n + c·Δt, so we want:
    // deriv(t_substage) = (1-c)·deriv(t^n) + c·deriv(t^{n+1})
    Real weight_n = 1.0 - c;
    Real weight_np1 = c;

    const Real* dx_L0 = geom_L0.CellSize();
    const Real* dx_L1 = geom_L1.CellSize();
    auto matches_two_to_one = [](Real coarse_dx, Real fine_dx) noexcept {
        Real expected = 2.0 * fine_dx;
        Real scale = std::max(std::abs(coarse_dx), std::abs(expected));
        return std::abs(coarse_dx - expected) <= 1.0e-12 * scale;
    };
    auto matches_one_to_one = [](Real coarse_dx, Real fine_dx) noexcept {
        Real scale = std::max(std::abs(coarse_dx), std::abs(fine_dx));
        return std::abs(coarse_dx - fine_dx) <= 1.0e-12 * scale;
    };

    if (!matches_two_to_one(dx_L0[0], dx_L1[0]) ||
        !matches_two_to_one(dx_L0[1], dx_L1[1]) ||
        !matches_one_to_one(dx_L0[2], dx_L1[2])) {
        amrex::Abort(
            "interpolateCoreDerivativeToL1Stage currently supports only ref_ratio = (2,2,1) between L0 and L1"
        );
    }

    // The prolongation kernels below are implemented only for ref_ratio = (2,2,1).
    IntVect ref_ratio(2, 2, 1);

    // Blend the coarse derivatives to the RK4 substage time first, then reuse
    // the same configurable spatial prolongation machinery as eta/phi.
    for (MFIter mfi(deriv_stage_scratch, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto deriv_stage_arr = deriv_stage_scratch.array(mfi);
        auto deriv_n_arr = deriv_L0_n.const_array(mfi);
        auto deriv_np1_arr = deriv_L0_np1.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            deriv_stage_arr(i,j,k) =
                weight_n * deriv_n_arr(i,j,k) + weight_np1 * deriv_np1_arr(i,j,k);
        });
    }

    crescent::interpolation::interpolateScalarParentToChild(
        deriv_stage_scratch,
        deriv_interp,
        geom_L0,
        geom_L1,
        ref_ratio,
        coarse_patch_scratch,
        interpolation_order
    );

    deriv_interp.FillBoundary(geom_L1.periodicity());
}

void blendL1DerivativeWithCore(
    const MultiFab& native_l1_derivative,
    const MultiFab& core_derivative_on_l1,
    const MultiFab& sigma,
    MultiFab& blended_derivative,
    MultiFab& correction_delta_raw,
    MultiFab& correction_delta_smooth,
    int correction_y_smoothing_passes)
{
    // Blend the native L1 derivative with the L0-on-L1 interpolation:
    // B[native_l1, core_on_l1] = (1-σ)·native_l1 + σ·core_on_l1
    //
    // Optional y-smoothing is applied only to the interpolation correction
    // (core_on_l1 - native_l1) inside the mixing zone, not to the full field.

    if (correction_y_smoothing_passes <= 0) {
        for (MFIter mfi(blended_derivative, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto native_arr = native_l1_derivative.const_array(mfi);
            auto core_arr = core_derivative_on_l1.const_array(mfi);
            auto sigma_arr = sigma.const_array(mfi);
            auto blended_arr = blended_derivative.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Real s = sigma_arr(i,j,k);
                blended_arr(i,j,k) = (1.0 - s) * native_arr(i,j,k) + s * core_arr(i,j,k);
            });
        }
        return;
    }

    for (MFIter mfi(correction_delta_raw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto correction_delta_arr = correction_delta_raw.array(mfi);
        auto native_arr = native_l1_derivative.const_array(mfi);
        auto core_arr = core_derivative_on_l1.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            correction_delta_arr(i,j,k) = core_arr(i,j,k) - native_arr(i,j,k);
        });
    }

    const Box global_box = blended_derivative.boxArray().minimalBox();
    const int j_lo = global_box.smallEnd(1);
    const int j_hi = global_box.bigEnd(1);

    MultiFab* current_correction_delta = &correction_delta_raw;
    MultiFab* next_correction_delta = &correction_delta_smooth;

    for (int pass = 0; pass < correction_y_smoothing_passes; ++pass) {
        current_correction_delta->FillBoundary();

        for (MFIter mfi(*next_correction_delta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto in_arr = current_correction_delta->const_array(mfi);
            auto out_arr = next_correction_delta->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const int jm = (j > j_lo) ? (j - 1) : j;
                const int jp = (j < j_hi) ? (j + 1) : j;
                out_arr(i,j,k) = 0.25 * in_arr(i,jm,k)
                               + 0.5 * in_arr(i,j,k)
                               + 0.25 * in_arr(i,jp,k);
            });
        }

        std::swap(current_correction_delta, next_correction_delta);
    }

    for (MFIter mfi(blended_derivative, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto native_arr = native_l1_derivative.const_array(mfi);
        auto core_arr = core_derivative_on_l1.const_array(mfi);
        auto sigma_arr = sigma.const_array(mfi);
        auto correction_delta_raw_arr = correction_delta_raw.const_array(mfi);
        auto correction_delta_smooth_arr = current_correction_delta->const_array(mfi);
        auto blended_arr = blended_derivative.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real s = sigma_arr(i,j,k);
            const Real base =
                (1.0 - s) * native_arr(i,j,k) + s * core_arr(i,j,k);
            const Real mixing_zone_correction =
                s * (1.0 - s)
                * (correction_delta_smooth_arr(i,j,k) - correction_delta_raw_arr(i,j,k));
            blended_arr(i,j,k) = base + mixing_zone_correction;
        });
    }
}

void computeSurfaceRHS_L1(
    MultiFab& rhs_eta,
    MultiFab& rhs_phi,
    const MultiFab& eta,
    const MultiFab& eta_x,
    const MultiFab& eta_y,
    const MultiFab& phi_x,
    const MultiFab& phi_y,
    const MultiFab& phi_z,
    int rk_stage,
    Real gravity)
{
    // Full nonlinear HOS equations (Zakharov-Craig-Sulem formulation),
    // consistent with the core-path computeSurfaceRHS in Source/core/CRESCENT_Dynamics.cpp:
    //
    //   ∂ₜη = -(∇ψ·∇η) + (1 + |∇η|²) · w
    //   ∂ₜψ = -g·η - ½|∇ψ|² + ½(1 + |∇η|²) · w²
    //
    // NOTE: No intermediate dealiasing here — L1 blended derivatives are NOT
    // periodic, so FFT-based dealiasing would introduce Gibbs artifacts.
    // Filtering is applied only to periodic data (C3-repaired copies in Step 1.5).
    //
    // CRITICAL: Write to component rk_stage (0-3) for proper RK4 update.
    for (MFIter mfi(rhs_eta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto eta_arr = eta.const_array(mfi);
        auto eta_x_arr = eta_x.const_array(mfi);
        auto eta_y_arr = eta_y.const_array(mfi);
        auto phi_x_arr = phi_x.const_array(mfi);
        auto phi_y_arr = phi_y.const_array(mfi);
        auto phi_z_arr = phi_z.const_array(mfi);
        auto rhs_eta_arr = rhs_eta.array(mfi);
        auto rhs_phi_arr = rhs_phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // metric = 1 + |∇η|²
            Real metric = 1.0 + eta_x_arr(i,j,k)*eta_x_arr(i,j,k)
                              + eta_y_arr(i,j,k)*eta_y_arr(i,j,k);

            // ∂ₜη = -(∇ψ·∇η) + metric · w
            rhs_eta_arr(i,j,k,rk_stage) =
                -(phi_x_arr(i,j,k) * eta_x_arr(i,j,k)
                + phi_y_arr(i,j,k) * eta_y_arr(i,j,k))
                + metric * phi_z_arr(i,j,k);

            // ∂ₜψ = -g·η - ½|∇ψ|² + ½·metric·w²
            Real grad_phi_sq = phi_x_arr(i,j,k) * phi_x_arr(i,j,k)
                             + phi_y_arr(i,j,k) * phi_y_arr(i,j,k);

            rhs_phi_arr(i,j,k,rk_stage) = -gravity * eta_arr(i,j,k)
                                        - 0.5 * grad_phi_sq
                                        + 0.5 * metric * phi_z_arr(i,j,k) * phi_z_arr(i,j,k);
        });
    }
}

} // namespace crescent::nested::rk4

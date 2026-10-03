#ifndef CRESCENT_NESTED_TYPES_H
#define CRESCENT_NESTED_TYPES_H

#include <algorithm>
#include <array>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_FreeSurfaceViews.H"

namespace crescent::nested {

struct InitializationValidationScratch
{
    amrex::MultiFab eta_core_interp;
    amrex::MultiFab phi_core_interp;
    amrex::MultiFab eta_x_core_interp;
    amrex::MultiFab eta_y_core_interp;
    amrex::MultiFab phi_x_core_interp;
    amrex::MultiFab phi_y_core_interp;
    amrex::MultiFab phi_core_interp_reference;
    amrex::MultiFab eta_prepared_temp;
    amrex::MultiFab phi_prepared_temp;
    amrex::MultiFab eta_prepared_x_temp;
    amrex::MultiFab eta_prepared_y_temp;
    amrex::MultiFab phi_prepared_x_temp;
    amrex::MultiFab phi_prepared_y_temp;
    amrex::MultiFab phi_x_blended_temp;
    amrex::MultiFab phi_y_blended_temp;
    amrex::MultiFab coarse_interp_patch;
};

struct BoundaryDiagnostic
{
    int i_lo = 0;
    int i_hi = 0;
    int j_mid = 0;
    amrex::Real x_left = 0.0;
    amrex::Real x_right = 0.0;
    amrex::Real y_mid = 0.0;
    amrex::Real phi_left = 0.0;
    amrex::Real phi_right = 0.0;
    amrex::Real phi_left_theory = 0.0;
    amrex::Real phi_right_theory = 0.0;
};

struct StageDerivativeViews {
    crescent::shared::DerivativeBundleView core_n;
    crescent::shared::DerivativeBundleView core_np1;
    crescent::shared::DerivativeBundleView core_derivatives_on_l1;
    crescent::shared::DerivativeBundleView blended_derivatives;
};

struct OutputDerivativeViews {
    crescent::shared::DerivativeBundleView core_np1;
    crescent::shared::DerivativeBundleView core_derivatives_on_l1;
    crescent::shared::DerivativeBundleView output_derivatives;
};

struct PeriodicDerivativeView {
    amrex::MultiFab& eta_periodic;
    amrex::MultiFab& phi_periodic;
    amrex::MultiFab& eta_linear_prefilter_periodic;
    amrex::MultiFab& phi_linear_prefilter_periodic;

    crescent::shared::FreeSurfaceStateView periodicState()
    {
        return {eta_periodic, phi_periodic};
    }

    crescent::shared::ConstFreeSurfaceStateView periodicState() const
    {
        return {eta_periodic, phi_periodic};
    }

    crescent::shared::FreeSurfaceStateView linearPrefilterState()
    {
        return {
            eta_linear_prefilter_periodic,
            phi_linear_prefilter_periodic
        };
    }

    crescent::shared::ConstFreeSurfaceStateView linearPrefilterState() const
    {
        return {
            eta_linear_prefilter_periodic,
            phi_linear_prefilter_periodic
        };
    }
};

struct PeriodicRealBasisCache {
    int npts = 0;
    int nharmonics = 0;
    std::vector<std::vector<amrex::Real>> basis;
};

struct C3RepairWorkspace {
    // Rank-0 gather scratch used by applyC3InterfaceSeamRepair. For runs that
    // do not regrid after setup, this single-box layout can be built once and
    // reused for every substage repair.
    amrex::BoxArray ba_global;
    amrex::DistributionMapping dm_global;
    amrex::MultiFab field_global;

    // Precomputed C3 interface basis matrices (deterministic from k, nx, ny,
    // bx, by, dx, dy). Rebuilt only when those parameters change.
    std::vector<std::vector<amrex::Real>> qx_basis;
    std::vector<std::vector<amrex::Real>> qy_basis;
    std::vector<std::vector<amrex::Real>> bx_basis;
    std::vector<std::vector<amrex::Real>> by_basis;
    PeriodicRealBasisCache x_pass_transverse_basis;
    PeriodicRealBasisCache x_pass_projection_basis;
    PeriodicRealBasisCache y_pass_transverse_basis;
    PeriodicRealBasisCache y_pass_projection_basis;
    int basis_k = 0;
    int basis_nx = 0;
    int basis_ny = 0;
    int basis_bx = 0;
    int basis_by = 0;
    amrex::Real basis_dx = 0.0;
    amrex::Real basis_dy = 0.0;

    // Per-pass line scratch. Sized once at setup from (nx, ny, k) and reused
    // on every substage repair call. x_pass buffers assume ny rows,
    // y_pass buffers assume nx rows. u_line is shared between passes and
    // sized to max(nx, ny).
    std::vector<amrex::Real> f_2d_scratch;
    std::vector<amrex::Real> x_correction;
    std::vector<amrex::Real> x_alpha;
    std::vector<amrex::Real> x_alpha_work;
    std::vector<std::array<amrex::Real, 4>> x_deltas;
    std::vector<amrex::Real> x_target_left;
    std::vector<amrex::Real> x_target_right;
    std::vector<amrex::Real> x_ref;
    std::vector<amrex::Real> x_u_left;
    std::vector<amrex::Real> x_u_right;
    std::vector<amrex::Real> x_correction_work;
    std::vector<amrex::Real> y_alpha;
    std::vector<std::array<amrex::Real, 4>> y_deltas;
    std::vector<amrex::Real> y_target_bottom;
    std::vector<amrex::Real> y_target_top;
    std::vector<amrex::Real> u_line_scratch;

    C3RepairWorkspace() = default;

    void define(const amrex::Box& global_box)
    {
        ba_global = amrex::BoxArray(global_box);
        dm_global = amrex::DistributionMapping(ba_global);
        field_global.define(ba_global, dm_global, 1, 0);
    }

    void resizeLineScratches(int k, int nx, int ny)
    {
        const int max_dim = std::max(nx, ny);
        f_2d_scratch.resize(nx * ny);
        x_correction.resize(nx * ny);
        x_alpha.resize(2 * k * ny);
        x_alpha_work.resize(2 * k * ny);
        x_deltas.resize(ny);
        x_target_left.resize(ny);
        x_target_right.resize(ny);
        x_ref.resize(ny);
        x_u_left.resize(ny);
        x_u_right.resize(ny);
        x_correction_work.resize(nx * ny);
        y_alpha.resize(nx * 2 * k);
        y_deltas.resize(nx);
        y_target_bottom.resize(nx);
        y_target_top.resize(nx);
        u_line_scratch.resize(max_dim);
    }

    [[nodiscard]] bool basisMatches(
        int k, int nx, int ny, int bx, int by,
        amrex::Real dx, amrex::Real dy) const noexcept
    {
        return basis_k == k && basis_nx == nx && basis_ny == ny
            && basis_bx == bx && basis_by == by
            && basis_dx == dx && basis_dy == dy;
    }
};

struct PostRK4StateScratch {
    // Persistent scratch used once per step by applyPostRK4CoreStateBlend to
    // hold (a) the filtered periodic L1 state (post C3-repair and filter) and
    // (b) the L0 state interpolated onto the L1 grid prior to sigma-blending.
    // The coarse interpolation patch reuses the same layout as the
    // interpolateParentToChild gather/prolong path so post-RK4 state blending
    // stops allocating a coarse patch every timestep.
    amrex::MultiFab periodic_l1_eta;
    amrex::MultiFab periodic_l1_phi;
    amrex::MultiFab interpolated_core_eta;
    amrex::MultiFab interpolated_core_phi;
    amrex::MultiFab coarse_interp_patch;

    PostRK4StateScratch() = default;

    void define(const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm,
                int nghost,
                const amrex::IntVect& ref_ratio)
    {
        periodic_l1_eta.define(ba, dm, 1, nghost);
        periodic_l1_phi.define(ba, dm, 1, nghost);
        interpolated_core_eta.define(ba, dm, 1, nghost);
        interpolated_core_phi.define(ba, dm, 1, nghost);
        coarse_interp_patch.define(amrex::coarsen(ba, ref_ratio), dm, 1, 2);
    }

    crescent::shared::FreeSurfaceStateView periodicL1View()
    {
        return {periodic_l1_eta, periodic_l1_phi};
    }

    crescent::shared::ConstFreeSurfaceStateView periodicL1ConstView() const
    {
        return {periodic_l1_eta, periodic_l1_phi};
    }

    crescent::shared::FreeSurfaceStateView interpolatedCoreView()
    {
        return {interpolated_core_eta, interpolated_core_phi};
    }

    crescent::shared::ConstFreeSurfaceStateView interpolatedCoreConstView() const
    {
        return {interpolated_core_eta, interpolated_core_phi};
    }

    amrex::MultiFab& coarseInterpPatch()
    {
        return coarse_interp_patch;
    }
};

struct OmegaAEdgeDiagnosticScratch {
    amrex::MultiFab eta_root;
    amrex::MultiFab phi_root;
    amrex::MultiFab eta_x_core_interp_root;
    amrex::MultiFab eta_y_core_interp_root;
    amrex::MultiFab phi_x_core_interp_root;
    amrex::MultiFab phi_y_core_interp_root;
    amrex::MultiFab eta_x_rk_root;
    amrex::MultiFab eta_y_rk_root;
    amrex::MultiFab phi_x_root;
    amrex::MultiFab phi_y_root;
    amrex::MultiFab phi_x_rk_root;
    amrex::MultiFab phi_y_rk_root;
    amrex::MultiFab phi_z_rk_root;
    amrex::MultiFab sigma_root;

    void define(const amrex::BoxArray& root_ba,
                const amrex::DistributionMapping& root_dm)
    {
        eta_root.define(root_ba, root_dm, 1, 0);
        phi_root.define(root_ba, root_dm, 1, 0);
        eta_x_core_interp_root.define(root_ba, root_dm, 1, 0);
        eta_y_core_interp_root.define(root_ba, root_dm, 1, 0);
        phi_x_core_interp_root.define(root_ba, root_dm, 1, 0);
        phi_y_core_interp_root.define(root_ba, root_dm, 1, 0);
        eta_x_rk_root.define(root_ba, root_dm, 1, 0);
        eta_y_rk_root.define(root_ba, root_dm, 1, 0);
        phi_x_root.define(root_ba, root_dm, 1, 0);
        phi_y_root.define(root_ba, root_dm, 1, 0);
        phi_x_rk_root.define(root_ba, root_dm, 1, 0);
        phi_y_rk_root.define(root_ba, root_dm, 1, 0);
        phi_z_rk_root.define(root_ba, root_dm, 1, 0);
        sigma_root.define(root_ba, root_dm, 1, 0);
    }
};

struct PhiStageDiagnosticScratch {
    amrex::MultiFab phi_stage_snapshot;
    amrex::MultiFab phi_raw_root;
    amrex::MultiFab phi_c3_root;
    amrex::MultiFab phi_filtered_root;
    amrex::MultiFab phi_y_raw;
    amrex::MultiFab phi_y_c3;
    amrex::MultiFab phi_y_filtered;
    amrex::MultiFab phi_y_raw_root;
    amrex::MultiFab phi_y_c3_root;
    amrex::MultiFab phi_y_filtered_root;

    void define(const amrex::BoxArray& level_ba,
                const amrex::DistributionMapping& level_dm,
                const amrex::BoxArray& root_ba,
                const amrex::DistributionMapping& root_dm)
    {
        phi_stage_snapshot.define(level_ba, level_dm, 1, 0);
        phi_raw_root.define(root_ba, root_dm, 1, 0);
        phi_c3_root.define(root_ba, root_dm, 1, 0);
        phi_filtered_root.define(root_ba, root_dm, 1, 0);
        phi_y_raw.define(level_ba, level_dm, 1, 0);
        phi_y_c3.define(level_ba, level_dm, 1, 0);
        phi_y_filtered.define(level_ba, level_dm, 1, 0);
        phi_y_raw_root.define(root_ba, root_dm, 1, 0);
        phi_y_c3_root.define(root_ba, root_dm, 1, 0);
        phi_y_filtered_root.define(root_ba, root_dm, 1, 0);
    }
};

struct RuntimeDiagnosticScratch {
    OmegaAEdgeDiagnosticScratch omega_a_edge;
    PhiStageDiagnosticScratch phi_stage;

    void define(const amrex::BoxArray& level_ba,
                const amrex::DistributionMapping& level_dm,
                const amrex::BoxArray& root_ba,
                const amrex::DistributionMapping& root_dm)
    {
        omega_a_edge.define(root_ba, root_dm);
        phi_stage.define(level_ba, level_dm, root_ba, root_dm);
    }
};

struct PeriodicDerivativeScratch {
    amrex::MultiFab eta_periodic;
    amrex::MultiFab phi_periodic;
    amrex::MultiFab eta_linear_prefilter_periodic;
    amrex::MultiFab phi_linear_prefilter_periodic;
    // Throwaway scratch used by cleanupNestedL1StagePhiZ as the second
    // (phi-slot) argument of the paired Hou-Li filter when filtering phi_z
    // alone. Defined once per level, always zero-ghost.
    amrex::MultiFab phi_z_filter_scratch;

    PeriodicDerivativeScratch() = default;

    PeriodicDerivativeScratch(const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm,
                              int ngrow = 0)
    {
        define(ba, dm, ngrow);
    }

    void define(const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm,
                int ngrow = 0)
    {
        eta_periodic.define(ba, dm, 1, ngrow);
        phi_periodic.define(ba, dm, 1, ngrow);
        eta_linear_prefilter_periodic.define(ba, dm, 1, ngrow);
        phi_linear_prefilter_periodic.define(ba, dm, 1, ngrow);
        phi_z_filter_scratch.define(ba, dm, 1, 0);
    }

    PeriodicDerivativeView view()
    {
        return {eta_periodic,
                phi_periodic,
                eta_linear_prefilter_periodic,
                phi_linear_prefilter_periodic};
    }

    crescent::shared::FreeSurfaceStateView periodicState()
    {
        return {eta_periodic, phi_periodic};
    }

    crescent::shared::ConstFreeSurfaceStateView periodicState() const
    {
        return {eta_periodic, phi_periodic};
    }

    crescent::shared::FreeSurfaceStateView linearPrefilterState()
    {
        return {
            eta_linear_prefilter_periodic,
            phi_linear_prefilter_periodic
        };
    }

    crescent::shared::ConstFreeSurfaceStateView linearPrefilterState() const
    {
        return {
            eta_linear_prefilter_periodic,
            phi_linear_prefilter_periodic
        };
    }

    void snapshotLinearPrefilterCopies()
    {
        amrex::MultiFab::Copy(eta_linear_prefilter_periodic, eta_periodic, 0, 0, 1, 0);
        amrex::MultiFab::Copy(phi_linear_prefilter_periodic, phi_periodic, 0, 0, 1, 0);
    }
};

using LinearDerivativeSources = crescent::shared::ConstFreeSurfaceStateView;

struct ColumnExtrema {
    amrex::Real min_value;
    amrex::Real max_value;
    int min_j;
    int max_j;
    int sample_count;
};

struct BoundaryColumnExtrema {
    ColumnExtrema left;
    ColumnExtrema right;
};

} // namespace crescent::nested

#endif // CRESCENT_NESTED_TYPES_H

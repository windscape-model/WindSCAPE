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

namespace crescent::nested {

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
    int basis_k = 0;
    int basis_nx = 0;
    int basis_ny = 0;
    int basis_bx = 0;
    int basis_by = 0;
    amrex::Real basis_dx = 0.0;
    amrex::Real basis_dy = 0.0;

    // Host line scratch, sized once at setup: the gathered field and one
    // line of length max(nx, ny).
    std::vector<amrex::Real> f_2d_scratch;
    std::vector<amrex::Real> u_line_scratch;

    C3RepairWorkspace() = default;

    void define(const amrex::Box& global_box)
    {
        ba_global = amrex::BoxArray(global_box);
        dm_global = amrex::DistributionMapping(ba_global);
        field_global.define(ba_global, dm_global, 1, 0);
    }

    void resizeLineScratches(int nx, int ny)
    {
        f_2d_scratch.resize(nx * ny);
        u_line_scratch.resize(std::max(nx, ny));
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
    // The parent state prolonged to the child for the post-RK4 state
    // relaxation, and the coarse patch of that prolongation.
    amrex::MultiFab interpolated_core_eta;
    amrex::MultiFab interpolated_core_phi;
    amrex::MultiFab coarse_interp_patch;

    void define(const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm,
                int nghost,
                const amrex::IntVect& ref_ratio)
    {
        interpolated_core_eta.define(ba, dm, 1, nghost);
        interpolated_core_phi.define(ba, dm, 1, nghost);
        coarse_interp_patch.define(amrex::coarsen(ba, ref_ratio), dm, 1, 2);
    }

    amrex::MultiFab& coarseInterpPatch()
    {
        return coarse_interp_patch;
    }
};

// Rank-0 (pinned) copies of the fields that the per-step nested diagnostics
// read, in the column order of omega_a_edge_diagnostic.csv. The oscillation
// diagnostic reuses one entry for its sequential gathers.
struct RuntimeDiagnosticScratch {
    static constexpr int nfields = 12;
    std::array<amrex::MultiFab, nfields> field_root;

    void define(const amrex::BoxArray& root_ba,
                const amrex::DistributionMapping& root_dm)
    {
        for (auto& root : field_root) {
            root.define(root_ba, root_dm, 1, 0,
                        amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        }
    }
};

// Periodic copies of the child state (C3-repaired and filtered), the only
// child fields that enter spectral operations in the total-field coupling.
struct PeriodicDerivativeScratch {
    amrex::MultiFab eta_periodic;
    amrex::MultiFab phi_periodic;

    // True from the end-of-step output preparation until the next stage:
    // the copies and the native derivative bundle then hold the preparation
    // of the current child state.
    bool holds_step_end_preparation = false;

    void define(const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm)
    {
        eta_periodic.define(ba, dm, 1, 0);
        phi_periodic.define(ba, dm, 1, 0);
    }
};

} // namespace crescent::nested

#endif // CRESCENT_NESTED_TYPES_H

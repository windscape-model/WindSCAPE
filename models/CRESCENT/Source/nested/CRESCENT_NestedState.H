#ifndef CRESCENT_NESTED_STATE_H
#define CRESCENT_NESTED_STATE_H

#include <array>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_NestedDomains.H"
#include "CRESCENT_NestedTypes.H"

namespace crescent::nested {

struct DerivativeFieldStorage
{
    amrex::Vector<amrex::MultiFab> eta_x;
    amrex::Vector<amrex::MultiFab> eta_y;
    amrex::Vector<amrex::MultiFab> phi_x;
    amrex::Vector<amrex::MultiFab> phi_y;
    amrex::Vector<amrex::MultiFab> phi_z;

    void resize(int nlevs_max)
    {
        eta_x.resize(nlevs_max);
        eta_y.resize(nlevs_max);
        phi_x.resize(nlevs_max);
        phi_y.resize(nlevs_max);
        phi_z.resize(nlevs_max);
    }

    void defineLevel(
        int lev,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        int nghost)
    {
        eta_x[lev].define(ba, dm, 1, nghost);
        eta_y[lev].define(ba, dm, 1, nghost);
        phi_x[lev].define(ba, dm, 1, nghost);
        phi_y[lev].define(ba, dm, 1, nghost);
        phi_z[lev].define(ba, dm, 1, nghost);
    }

    crescent::shared::DerivativeBundleView bundleAt(int lev)
    {
        return {eta_x[lev], eta_y[lev], phi_x[lev], phi_y[lev], phi_z[lev]};
    }
};

// The parent derivative bundle at the end of the current step, the parent
// input of the child output bundles.
struct CoreDerivativeSnapshots
{
    DerivativeFieldStorage np1;

    void resize(int nlevs_max)
    {
        np1.resize(nlevs_max);
    }

    void defineLevel(
        int lev,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        int nghost)
    {
        np1.defineLevel(lev, ba, dm, nghost);
    }
};

struct CouplingZoneState
{
    amrex::Vector<domains::DomainStructure> domains;
    amrex::Vector<int> c3_buffer_width_x_computed;
    amrex::Vector<int> c3_buffer_width_y_computed;
    // Retained so current Tecplot/CRESCENT diagnostics can report the expected
    // mask field. The production nested path no longer populates damping masks.
    amrex::Vector<amrex::MultiFab> chi_D_mask;
    amrex::Vector<amrex::MultiFab> sigma_blend;

    void resize(int nlevs_max)
    {
        domains.resize(nlevs_max);
        c3_buffer_width_x_computed.resize(nlevs_max, 0);
        c3_buffer_width_y_computed.resize(nlevs_max, 0);
        chi_D_mask.resize(nlevs_max);
        sigma_blend.resize(nlevs_max);
    }

    void defineLevel(
        int lev,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm)
    {
        chi_D_mask[lev].define(ba, dm, 1, 0);
        sigma_blend[lev].define(ba, dm, 1, 0);
        chi_D_mask[lev].setVal(0.0);
    }
};

struct CouplingDerivativeState
{
    CoreDerivativeSnapshots core_snapshots;
    DerivativeFieldStorage core_derivatives_on_l1;
    DerivativeFieldStorage l1_blended_derivatives;
    // Parent-grid and coarse-patch scratch reused by the parent-to-child
    // prolongations.
    amrex::MultiFab core_stage_derivative_scratch;
    amrex::MultiFab core_interp_patch_scratch;
    amrex::MultiFab adjustment_w_lin_scratch; // w_lin workspace for the Dommermuth ramp on the child
    // Tendency-blend coupling (l1_coupling_mode == 1): per-stage native L1
    // RHS from the periodized fields and the interpolated L0 stage RHS.
    amrex::MultiFab rhs_eta_native;
    amrex::MultiFab rhs_phi_native;
    amrex::MultiFab rhs_eta_core_on_l1;
    amrex::MultiFab rhs_phi_core_on_l1;

    void resize(int nlevs_max)
    {
        core_snapshots.resize(nlevs_max);
        core_derivatives_on_l1.resize(nlevs_max);
        l1_blended_derivatives.resize(nlevs_max);
    }

    void defineLevel(
        int lev,
        const amrex::BoxArray& core_ba,
        const amrex::DistributionMapping& core_dm,
        const amrex::BoxArray& nested_ba,
        const amrex::DistributionMapping& nested_dm,
        const amrex::IntVect& ref_ratio)
    {
        core_snapshots.defineLevel(0, core_ba, core_dm, 1);
        core_derivatives_on_l1.defineLevel(lev, nested_ba, nested_dm, 0);
        l1_blended_derivatives.defineLevel(lev, nested_ba, nested_dm, 0);
        core_stage_derivative_scratch.define(core_ba, core_dm, 1, 0);
        core_interp_patch_scratch.define(
            amrex::coarsen(nested_ba, ref_ratio),
            nested_dm,
            1,
            2);
        rhs_eta_native.define(nested_ba, nested_dm, 1, 0);
        rhs_phi_native.define(nested_ba, nested_dm, 1, 0);
        adjustment_w_lin_scratch.define(nested_ba, nested_dm, 1, 0);
        rhs_eta_core_on_l1.define(nested_ba, nested_dm, 1, 0);
        rhs_phi_core_on_l1.define(nested_ba, nested_dm, 1, 0);
    }
};

// Correction-field coupling (l1_coupling_mode == 2): per-RK-stage snapshots
// of the parent state/derivative bundle plus two parent-grid spectral
// quantities, captured observationally during the parent advance
// (docs/R2_CORRECTION_FIELD_SPEC.md §4). Defined only in mode 2.
struct CoreStageSnapshotFields
{
    amrex::MultiFab eta;
    amrex::MultiFab phi;
    amrex::MultiFab eta_x;
    amrex::MultiFab eta_y;
    amrex::MultiFab phi_x;
    amrex::MultiFab phi_y;
    amrex::MultiFab phi_z;
    amrex::MultiFab w_lin;    // IFFT(|k| phi_hat) on the parent grid
    amrex::MultiFab lap_phi;  // spectral laplacian of phi (coefficient-projected when requested) on the parent grid
    amrex::MultiFab m0_d;     // parent-dealiased eta_x^2 + eta_y^2
    amrex::MultiFab w0sq_d;   // parent-dealiased phi_z^2

    void define(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    {
        eta.define(ba, dm, 1, 0);
        phi.define(ba, dm, 1, 0);
        eta_x.define(ba, dm, 1, 0);
        eta_y.define(ba, dm, 1, 0);
        phi_x.define(ba, dm, 1, 0);
        phi_y.define(ba, dm, 1, 0);
        phi_z.define(ba, dm, 1, 0);
        w_lin.define(ba, dm, 1, 0);
        lap_phi.define(ba, dm, 1, 0);
        m0_d.define(ba, dm, 1, 0);
        w0sq_d.define(ba, dm, 1, 0);
    }
};

// Child-grid scratch for the correction-field stage: interpolated background
// fields, the tapered correction pair, and term scratch. Defined only in
// mode 2.
struct CorrectionFieldScratch
{
    static constexpr int num_stages = 4;
    std::array<CoreStageSnapshotFields, num_stages> core_stage;
    // Final parent background of the step (state u_L0^{n+1} and its
    // derivative bundle), captured for the output-bundle recompute so the
    // diagnostic phi_z carries the same order-2 background-correction cross
    // terms as the stage tendency (added 2026-08-16).
    CoreStageSnapshotFields core_final;
    bool core_stage_defined = false;

    // Background interpolated to the child at the current stage.
    CoreStageSnapshotFields bg_on_l1;
    // Tapered correction pair (chi * delta), FFT-admissible by compactness.
    amrex::MultiFab de;
    amrex::MultiFab dp;
    // Correction vertical velocity, metric increment, and term scratch.
    amrex::MultiFab dw;
    amrex::MultiFab dm_raw;
    amrex::MultiFab dm;
    amrex::MultiFab term_a;
    amrex::MultiFab term_b;
    amrex::MultiFab cross_a;
    amrex::MultiFab cross_b;

    void defineCoreStages(
        const amrex::BoxArray& core_ba,
        const amrex::DistributionMapping& core_dm)
    {
        for (auto& s : core_stage) {
            s.define(core_ba, core_dm);
        }
        core_final.define(core_ba, core_dm);
        core_stage_defined = true;
    }

    void defineChild(
        const amrex::BoxArray& nested_ba,
        const amrex::DistributionMapping& nested_dm)
    {
        bg_on_l1.define(nested_ba, nested_dm);
        de.define(nested_ba, nested_dm, 1, 0);
        dp.define(nested_ba, nested_dm, 1, 0);
        dw.define(nested_ba, nested_dm, 1, 0);
        dm_raw.define(nested_ba, nested_dm, 1, 0);
        dm.define(nested_ba, nested_dm, 1, 0);
        term_a.define(nested_ba, nested_dm, 1, 0);
        term_b.define(nested_ba, nested_dm, 1, 0);
        cross_a.define(nested_ba, nested_dm, 1, 0);
        cross_b.define(nested_ba, nested_dm, 1, 0);
    }
};

struct RuntimeState
{
    CouplingZoneState coupling_zones;
    CouplingDerivativeState coupling_derivatives;
    CorrectionFieldScratch correction_field;
    amrex::Vector<PeriodicDerivativeScratch> periodic_derivative_scratch;
    amrex::Vector<PostRK4StateScratch> post_rk4_state_scratch;
    amrex::Vector<C3RepairWorkspace> c3_repair_workspace;
    amrex::Vector<RuntimeDiagnosticScratch> runtime_diagnostic_scratch;

    void resize(int nlevs_max)
    {
        coupling_zones.resize(nlevs_max);
        coupling_derivatives.resize(nlevs_max);
        periodic_derivative_scratch.resize(nlevs_max);
        post_rk4_state_scratch.resize(nlevs_max);
        c3_repair_workspace.resize(nlevs_max);
        runtime_diagnostic_scratch.resize(nlevs_max);
    }
};

} // namespace crescent::nested

#endif // CRESCENT_NESTED_STATE_H

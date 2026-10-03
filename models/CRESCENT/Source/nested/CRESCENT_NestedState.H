#ifndef CRESCENT_NESTED_STATE_H
#define CRESCENT_NESTED_STATE_H

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

struct CoreDerivativeSnapshots
{
    DerivativeFieldStorage n;
    DerivativeFieldStorage np1;

    void resize(int nlevs_max)
    {
        n.resize(nlevs_max);
        np1.resize(nlevs_max);
    }

    void defineLevel(
        int lev,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        int nghost)
    {
        n.defineLevel(lev, ba, dm, nghost);
        np1.defineLevel(lev, ba, dm, nghost);
    }
};

struct CouplingZoneState
{
    amrex::Vector<domains::DomainStructure> domains;
    amrex::Vector<int> c3_buffer_width_computed;
    // Retained so current Tecplot/CRESCENT diagnostics can report the expected
    // mask field. The production nested path no longer populates damping masks.
    amrex::Vector<amrex::MultiFab> chi_D_mask;
    amrex::Vector<amrex::MultiFab> sigma_blend;

    void resize(int nlevs_max)
    {
        domains.resize(nlevs_max);
        c3_buffer_width_computed.resize(nlevs_max, 0);
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
    // Reused sequentially while blending L0 snapshots to an RK stage time and
    // prolonging the result onto the L1 grid. This keeps the active derivative
    // coupling path from reallocating the same core-grid and coarse-patch
    // scratch on every component interpolation.
    amrex::MultiFab core_stage_derivative_scratch;
    amrex::MultiFab core_interp_patch_scratch;
    // Reused sequentially by optional y-smoothing of the interpolation
    // correction during derivative blending.
    amrex::MultiFab blend_correction_delta_raw;
    amrex::MultiFab blend_correction_delta_smooth;

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
        const amrex::IntVect& ref_ratio,
        int correction_y_smoothing_passes)
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
        if (correction_y_smoothing_passes > 0) {
            blend_correction_delta_raw.define(nested_ba, nested_dm, 1, 1);
            blend_correction_delta_smooth.define(nested_ba, nested_dm, 1, 1);
        }
    }
};

struct RuntimeState
{
    CouplingZoneState coupling_zones;
    CouplingDerivativeState coupling_derivatives;
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

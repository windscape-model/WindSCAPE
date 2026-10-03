#ifndef CRESCENT_NESTED_PERIODIC_RESEARCH_DIAGNOSTICS_H
#define CRESCENT_NESTED_PERIODIC_RESEARCH_DIAGNOSTICS_H

#include <algorithm>
#include <cmath>
#include <fstream>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_NestedResearchDiagnostics.H"

namespace crescent::nested::research_diagnostics {

using LevelMeshView = crescent::shared::LevelMeshView;
using PeriodicDerivativeScratch = crescent::nested::PeriodicDerivativeScratch;

template <typename ApplyNestedC3RepairFn>
inline void runInitialC3RepairEffectDiagnostic(
    const LevelMeshView& level_mesh,
    const amrex::MultiFab& eta,
    const amrex::MultiFab& phi,
    PeriodicDerivativeScratch& periodic_fields,
    ApplyNestedC3RepairFn&& apply_nested_c3_repair_pair)
{
    amrex::MultiFab::Copy(periodic_fields.eta_periodic, eta, 0, 0, 1, 0);
    amrex::MultiFab::Copy(periodic_fields.phi_periodic, phi, 0, 0, 1, 0);

    const amrex::Real max_dy_before =
        computeMaxAdjacentYDifference(periodic_fields.eta_periodic);

    amrex::MultiFab phi_before_c3(level_mesh.grids, level_mesh.dmap, 1, 0);
    amrex::MultiFab::Copy(
        phi_before_c3,
        periodic_fields.phi_periodic,
        0,
        0,
        1,
        0);
    const amrex::Real phi_before_min = phi_before_c3.min(0);
    const amrex::Real phi_before_max = phi_before_c3.max(0);

    apply_nested_c3_repair_pair(
        periodic_fields.eta_periodic,
        periodic_fields.phi_periodic,
        0);

    const amrex::Real max_dy_after =
        computeMaxAdjacentYDifference(periodic_fields.eta_periodic);
    const amrex::Real phi_after_min = periodic_fields.phi_periodic.min(0);
    const amrex::Real phi_after_max = periodic_fields.phi_periodic.max(0);

    amrex::MultiFab phi_correction(level_mesh.grids, level_mesh.dmap, 1, 0);
    amrex::MultiFab::Copy(
        phi_correction,
        periodic_fields.phi_periodic,
        0,
        0,
        1,
        0);
    amrex::MultiFab::Subtract(phi_correction, phi_before_c3, 0, 0, 1, 0);
    const amrex::Real corr_min = phi_correction.min(0);
    const amrex::Real corr_max = phi_correction.max(0);

    amrex::Print() << "\n=== C³-Repair Effect Diagnostic ===" << std::endl;
    amrex::Print() << "  Eta max |Δy| before C³-repair: " << max_dy_before
                   << std::endl;
    amrex::Print() << "  Eta max |Δy| after C³-repair:  " << max_dy_after
                   << std::endl;
    amrex::Print() << "  Phi before C³-repair: [" << phi_before_min << ", "
                   << phi_before_max << "]" << std::endl;
    amrex::Print() << "  Phi after C³-repair:  [" << phi_after_min << ", "
                   << phi_after_max << "]" << std::endl;
    amrex::Print() << "  Correction applied:   [" << corr_min << ", "
                   << corr_max << "]" << std::endl;
    amrex::Print() << "  Max correction magnitude: "
                   << std::max(std::abs(corr_min), std::abs(corr_max))
                   << std::endl;

    std::ofstream csv_file("c3_repair_diagnostic.csv");
    csv_file << "i,j,x,y,phi_before,phi_after,correction\n";

    for (amrex::MFIter mfi(phi_before_c3); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const& phi_before_arr = phi_before_c3.array(mfi);
        auto const& phi_after_arr = periodic_fields.phi_periodic.array(mfi);
        auto const& phi_correction_arr = phi_correction.array(mfi);

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                const amrex::Real x =
                    level_mesh.geom.ProbLo(0) +
                    (i + 0.5) * level_mesh.geom.CellSize(0);
                const amrex::Real y =
                    level_mesh.geom.ProbLo(1) +
                    (j + 0.5) * level_mesh.geom.CellSize(1);
                csv_file << i << "," << j << "," << x << "," << y << ","
                         << phi_before_arr(i, j, 0) << ","
                         << phi_after_arr(i, j, 0) << ","
                         << phi_correction_arr(i, j, 0) << "\n";
            }
        }
    }

    amrex::Print() << "  Saved detailed data to c3_repair_diagnostic.csv"
                   << std::endl;
}

inline void logNestedL1PhiZCleanupDiagnostics(const amrex::MultiFab& eta_x,
                                              const amrex::MultiFab& phi_z)
{
    const amrex::Real eta_x_amp = eta_x.norm0();
    const amrex::Real phi_z_amp = phi_z.norm0();

    amrex::Print() << "  Step 2.5 (phi_z periodic cleanup): phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0)
                   << ", |eta_x|_inf=" << eta_x_amp
                   << ", |phi_z|_inf=" << phi_z_amp << std::endl;
}

inline void logNestedPreRepairYVariation(const amrex::MultiFab& eta_level,
                                         int verbosity)
{
    if (verbosity <= 1) {
        return;
    }

    amrex::Print() << "    eta[lev] max |Δy|: "
                   << computeMaxAdjacentYDifference(eta_level)
                   << " (before C³-repair)" << std::endl;
}

inline void logNestedPostRepairYVariation(const amrex::MultiFab& eta_periodic,
                                          int verbosity)
{
    if (verbosity <= 1) {
        return;
    }

    amrex::Print() << "    eta_periodic max |Δy|: "
                   << computeMaxAdjacentYDifference(eta_periodic)
                   << " (after C³-repair)" << std::endl;
}

inline void logNestedPeriodicPrefilterState(
    const amrex::MultiFab& eta_periodic,
    int log_stage,
    int verbosity)
{
    if (!(verbosity >= 2 && log_stage == 0)) {
        return;
    }

    amrex::Print() << "  Dealiasing periodic copies before derivative computation..."
                   << std::endl;
    amrex::Print() << "    eta_periodic BEFORE filter: min="
                   << eta_periodic.min(0)
                   << ", max=" << eta_periodic.max(0) << std::endl;
}

inline void logNestedPeriodicPostfilterState(
    const amrex::MultiFab& eta_periodic,
    int log_stage,
    int verbosity)
{
    if (!(verbosity >= 2 && log_stage == 0)) {
        return;
    }

    amrex::Print() << "    eta_periodic AFTER filter: min="
                   << eta_periodic.min(0)
                   << ", max=" << eta_periodic.max(0) << std::endl;
}

} // namespace crescent::nested::research_diagnostics

#endif // CRESCENT_NESTED_PERIODIC_RESEARCH_DIAGNOSTICS_H

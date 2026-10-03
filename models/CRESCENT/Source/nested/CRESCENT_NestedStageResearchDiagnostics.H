#ifndef CRESCENT_NESTED_STAGE_RESEARCH_DIAGNOSTICS_H
#define CRESCENT_NESTED_STAGE_RESEARCH_DIAGNOSTICS_H

#include <cmath>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_NestedResearchDiagnostics.H"

namespace crescent::nested::research_diagnostics {

inline bool shouldLogDetailedStageDiagnostics(int verbosity)
{
    return verbosity >= 2;
}

inline int centerColumnIndex(const amrex::MultiFab& field)
{
    const amrex::Box domain = field.boxArray().minimalBox();
    return (domain.smallEnd(0) + domain.bigEnd(0)) / 2;
}

inline amrex::Real cellCenterX(const amrex::Geometry& geom, int i)
{
    return geom.ProbLo(0) +
           (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * geom.CellSize(0);
}

inline void logNestedL1StageEntryDiagnostics(const amrex::MultiFab& eta_level,
                                             const amrex::Box& accurate_region,
                                             int irk,
                                             int verbosity)
{
    if (!shouldLogDetailedStageDiagnostics(verbosity)) {
        return;
    }

    const amrex::Real eta_amp_before = eta_level.norm0();
    const amrex::Real eta_rms_before =
        eta_level.norm2() / std::sqrt(eta_level.boxArray().numPts());

    amrex::Print() << "  Before stage " << irk << ": |eta|_inf=" << eta_amp_before
                   << ", |eta|_rms=" << eta_rms_before << std::endl;

    const auto boundary_extrema =
        sampleBoundaryColumnExtrema(
            eta_level,
            accurate_region.smallEnd(0),
            accurate_region.bigEnd(0),
            0);
    const auto& left_extrema = boundary_extrema.left;
    const auto& right_extrema = boundary_extrema.right;

    amrex::Print() << "    eta[lev] boundary y-variation (BEFORE operations): left Δ="
                   << columnVariation(left_extrema)
                   << ", right Δ=" << columnVariation(right_extrema) << std::endl;

    if (irk == 0 && left_extrema.sample_count > 0) {
        amrex::Print() << "      Left Omega_A boundary (i="
                       << accurate_region.smallEnd(0) << "): min="
                       << left_extrema.min_value
                       << ", max=" << left_extrema.max_value
                       << ", count=" << left_extrema.sample_count << std::endl;
        amrex::Print() << "      Right Omega_A boundary (i="
                       << accurate_region.bigEnd(0) << "): min="
                       << right_extrema.min_value
                       << ", max=" << right_extrema.max_value
                       << ", count=" << right_extrema.sample_count << std::endl;
    }
}

inline void logNestedL1CorePhiZReferenceDiagnostics(const amrex::Geometry& core_geom,
                                                    const amrex::MultiFab& core_phi_z,
                                                    const amrex::Geometry& level_geom,
                                                    const amrex::Box& accurate_region)
{
    const amrex::Real* prob_lo_core = core_geom.ProbLo();
    const amrex::Real* dx_arr_core = core_geom.CellSize();
    const int level_center_i = centerColumnIndex(core_phi_z);
    const amrex::Real center_x = cellCenterX(core_geom, level_center_i);
    const int i_sample_core = level_center_i;
    const auto sample_extrema =
        sampleColumnExtrema(core_phi_z, i_sample_core, 0);

    amrex::Print() << "    L0 phi_z center-column y-variation (i="
                   << i_sample_core << ", x=" << center_x << "): Δ="
                   << columnVariation(sample_extrema) << std::endl;

    const amrex::Real nested_left_boundary_x =
        cellCenterX(level_geom, accurate_region.smallEnd(0));
    const amrex::Real nested_right_boundary_x =
        cellCenterX(level_geom, accurate_region.bigEnd(0));

    const int core_left_boundary_i =
        static_cast<int>((nested_left_boundary_x - prob_lo_core[0]) / dx_arr_core[0]);
    const int core_right_boundary_i =
        static_cast<int>((nested_right_boundary_x - prob_lo_core[0]) / dx_arr_core[0]);
    const auto boundary_extrema =
        sampleBoundaryColumnExtrema(core_phi_z, core_left_boundary_i, core_right_boundary_i, 0);
    const auto& left_extrema = boundary_extrema.left;
    const auto& right_extrema = boundary_extrema.right;

    amrex::Print() << "    L0 phi_z at Omega_A x-boundaries: left (i="
                   << core_left_boundary_i << ", x=" << nested_left_boundary_x
                   << ") Δy="
                   << columnVariation(left_extrema)
                   << ", right (i=" << core_right_boundary_i
                   << ", x=" << nested_right_boundary_x << ") Δy="
                   << columnVariation(right_extrema) << std::endl;
}

inline void logNestedL1InterpolatedDerivativeDiagnostics(
    crescent::shared::DerivativeBundleView core_derivatives_on_l1,
    const amrex::Box& accurate_region)
{
    const amrex::Real eta_x_core_interp_amp = core_derivatives_on_l1.eta_x.norm0();
    const amrex::Real eta_y_core_interp_amp = core_derivatives_on_l1.eta_y.norm0();
    const amrex::Real phi_z_core_interp_amp = core_derivatives_on_l1.phi_z.norm0();

    amrex::Print() << "  Step 3 (L0 on L1): |eta_x_interp|_inf=" << eta_x_core_interp_amp
                   << ", |eta_y_interp|_inf=" << eta_y_core_interp_amp
                   << ", |phi_z_interp|_inf=" << phi_z_core_interp_amp << std::endl;

    const auto boundary_extrema =
        sampleBoundaryColumnExtrema(
            core_derivatives_on_l1.phi_z,
            accurate_region.smallEnd(0),
            accurate_region.bigEnd(0),
            0);
    const auto& left_extrema = boundary_extrema.left;
    const auto& right_extrema = boundary_extrema.right;

    amrex::Print() << "    phi_z_interp (L0 on L1) boundary y-variation: left Δ="
                   << columnVariation(left_extrema)
                   << " (min@j=" << left_extrema.min_j
                   << ", max@j=" << left_extrema.max_j << ")"
                   << ", right Δ=" << columnVariation(right_extrema)
                   << " (min@j=" << right_extrema.min_j
                   << ", max@j=" << right_extrema.max_j << ")" << std::endl;
}

inline void logNestedL1BlendedDerivativeDiagnostics(
    crescent::shared::DerivativeBundleView blended_bundle)
{
    const amrex::Real eta_x_RK_amp = blended_bundle.eta_x.norm0();
    const amrex::Real eta_y_RK_amp = blended_bundle.eta_y.norm0();
    const amrex::Real phi_z_RK_amp = blended_bundle.phi_z.norm0();

    amrex::Print() << "  Step 4 (Blend derivatives): |eta_x_RK|_inf=" << eta_x_RK_amp
                   << ", |eta_y_RK|_inf=" << eta_y_RK_amp
                   << ", |phi_z_RK|_inf=" << phi_z_RK_amp << std::endl;
}

inline void logNestedL1RhsDiagnostics(const amrex::MultiFab& rhs_eta_level,
                                      const amrex::MultiFab& rhs_phi_level,
                                      const amrex::Box& accurate_region,
                                      int irk)
{
    const amrex::Real rhs_eta_amp = rhs_eta_level.norm0();
    const amrex::Real rhs_phi_amp = rhs_phi_level.norm0();

    amrex::Print() << "  Step 5 (RHS from blended derivatives): |rhs_eta|_inf=" << rhs_eta_amp
                   << ", |rhs_phi|_inf=" << rhs_phi_amp << std::endl;

    const auto boundary_extrema =
        sampleBoundaryColumnExtrema(
            rhs_eta_level,
            accurate_region.smallEnd(0),
            accurate_region.bigEnd(0),
            irk);
    const auto& left_extrema = boundary_extrema.left;
    const auto& right_extrema = boundary_extrema.right;

    amrex::Print() << "    rhs_eta boundary y-variation: left Δ="
                   << columnVariation(left_extrema)
                   << ", right Δ=" << columnVariation(right_extrema)
                   << std::endl;
}

} // namespace crescent::nested::research_diagnostics

#endif // CRESCENT_NESTED_STAGE_RESEARCH_DIAGNOSTICS_H

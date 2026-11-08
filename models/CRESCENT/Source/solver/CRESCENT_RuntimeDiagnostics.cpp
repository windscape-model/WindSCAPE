#include "solver/CRESCENT_RuntimeDiagnostics.H"

// Per-step nested diagnostics: the Omega_A edge probes (a regression-gate
// file) and the oscillation and periodicity-seam metrics.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

#include <AMReX_FileSystem.H>

#include "CRESCENT_SolverContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_Utils.H"

using namespace amrex;

namespace {

std::string runtimeOutputDirectory(const std::string& plot_file)
{
    const std::size_t slash = plot_file.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string(".") : plot_file.substr(0, slash);
}

} // namespace

namespace crescent::solver::runtime_diagnostics {

void writeOmegaAEdgeDiagnostics(
    crescent::shared::LevelContext& level_context,
    amrex::Real time_next)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    auto& runtime = level_context.runtime();
    auto& reporting = runtime.reporting;
    const auto& coupling_zones = level_context.nestedRuntime().coupling_zones;
    const auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;
    auto& diagnostic_scratch =
        level_context.nestedRuntime().runtime_diagnostic_scratch[lev];
    const auto& level_geom = level_context.levelGeometry();

    if (!crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev) ||
        reporting.write_omega_a_edge_diagnostics != 1) {
        return;
    }
    if (runtime.istep[lev] % reporting.omega_a_edge_diagnostic_int != 0) {
        return;
    }

    const auto& domain = coupling_zones.domains[lev];
    const Box& box_A = domain.box_A;
    const Box& level_domain = level_geom.Domain();
    const int j_lo = box_A.smallEnd(1);
    const int j_hi = box_A.bigEnd(1);

    if (j_lo > j_hi) {
        return;
    }

    const int i_candidates[] = {
        box_A.smallEnd(0) - 1,
        box_A.smallEnd(0),
        box_A.smallEnd(0) + 1,
        box_A.bigEnd(0) - 1,
        box_A.bigEnd(0),
        box_A.bigEnd(0) + 1
    };
    const char* labels[] = {
        "left_outer",
        "left_edge",
        "left_inner",
        "right_inner",
        "right_edge",
        "right_outer"
    };

    // The probed fields, in the column order of the file: sigma, the child
    // state, the interpolated parent derivatives, and the blended child
    // derivatives.
    const auto& parent_on_child = coupling_derivatives.core_derivatives_on_l1;
    const auto& blended = coupling_derivatives.l1_blended_derivatives;
    const MultiFab* const fields[crescent::nested::RuntimeDiagnosticScratch::nfields] = {
        &coupling_zones.sigma_blend[lev],
        &solution.eta[lev],
        &solution.phi[lev],
        &parent_on_child.eta_x[lev],
        &parent_on_child.eta_y[lev],
        &parent_on_child.phi_x[lev],
        &parent_on_child.phi_y[lev],
        &blended.eta_x[lev],
        &blended.eta_y[lev],
        &blended.phi_x[lev],
        &blended.phi_y[lev],
        &blended.phi_z[lev]
    };
    constexpr int nfields = crescent::nested::RuntimeDiagnosticScratch::nfields;
    for (int f = 0; f < nfields; ++f) {
        crescent::utils::gatherToRoot(*fields[f], diagnostic_scratch.field_root[f]);
    }

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string out_dir = runtimeOutputDirectory(reporting.plot_file);
    const std::string filename = out_dir + "/omega_a_edge_diagnostic.csv";

    const bool initialize_main =
        !reporting.omega_a_edge_diagnostics_initialized ||
        !amrex::FileSystem::Exists(filename);
    std::ofstream ofs(
        filename,
        initialize_main ? std::ios::trunc : std::ios::app);
    if (!ofs.is_open()) {
        amrex::Abort("Failed to open omega_a_edge_diagnostic.csv for writing");
    }
    if (initialize_main) {
        ofs << "step,time,probe_label,i_probe,x_probe,j_lo,j_hi,"
               "sigma_min,sigma_max,eta_min,eta_max,phi_min,phi_max,"
               "eta_x_interp_min,eta_x_interp_max,eta_y_interp_min,eta_y_interp_max,"
               "phi_x_interp_min,phi_x_interp_max,phi_y_interp_min,phi_y_interp_max,"
               "eta_x_RK_min,eta_x_RK_max,eta_y_RK_min,eta_y_RK_max,"
               "phi_x_RK_min,phi_x_RK_max,phi_y_RK_min,phi_y_RK_max,"
               "phi_z_RK_min,phi_z_RK_max\n";
    }
    reporting.omega_a_edge_diagnostics_initialized = true;

    const Real* dx = level_geom.CellSize();
    const auto prob_lo = level_geom.ProbLoArray();

    for (int probe = 0; probe < 6; ++probe) {
        const int i_probe = i_candidates[probe];
        if (i_probe < level_domain.smallEnd(0) || i_probe > level_domain.bigEnd(0)) {
            continue;
        }

        const Real x_probe =
            prob_lo[0] + (static_cast<Real>(i_probe) + Real(0.5)) * dx[0];
        ofs << runtime.istep[lev] << "," << time_next << ","
            << labels[probe] << "," << i_probe << "," << x_probe << ","
            << j_lo << "," << j_hi;
        for (int f = 0; f < nfields; ++f) {
            const auto arr = diagnostic_scratch.field_root[f].const_arrays()[0];
            Real field_min = arr(i_probe, j_lo, 0);
            Real field_max = field_min;
            for (int j = j_lo + 1; j <= j_hi; ++j) {
                const Real val = arr(i_probe, j, 0);
                field_min = std::min(field_min, val);
                field_max = std::max(field_max, val);
            }
            ofs << "," << field_min << "," << field_max;
        }
        ofs << "\n";
    }
}

namespace {

struct SeamMetrics
{
    Real step_jump = 0.0;   // max |u(first) - u(last)| across the wrap
    Real d2_jump = 0.0;     // max |u(first+1) - 2 u(first) + u(last)| across the wrap
};

// Periodicity-seam mismatch of a root-gathered field. For an exactly periodic
// grid function, column 0 is the wrap-neighbor of column n-1: the step jump
// measures C0 mismatch, the wrap second difference measures slope mismatch.
SeamMetrics seamMetricsX(const amrex::Array4<const Real>& arr, const Box& bx)
{
    SeamMetrics m;
    const int i_lo = bx.smallEnd(0);
    const int i_hi = bx.bigEnd(0);
    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        const Real wrap_step = arr(i_lo, j, 0) - arr(i_hi, j, 0);
        const Real wrap_d2 =
            arr(i_lo + 1, j, 0) - 2.0 * arr(i_lo, j, 0) + arr(i_hi, j, 0);
        m.step_jump = std::max(m.step_jump, std::abs(wrap_step));
        m.d2_jump = std::max(m.d2_jump, std::abs(wrap_d2));
    }
    return m;
}

SeamMetrics seamMetricsY(const amrex::Array4<const Real>& arr, const Box& bx)
{
    SeamMetrics m;
    const int j_lo = bx.smallEnd(1);
    const int j_hi = bx.bigEnd(1);
    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
        const Real wrap_step = arr(i, j_lo, 0) - arr(i, j_hi, 0);
        const Real wrap_d2 =
            arr(i, j_lo + 1, 0) - 2.0 * arr(i, j_lo, 0) + arr(i, j_hi, 0);
        m.step_jump = std::max(m.step_jump, std::abs(wrap_step));
        m.d2_jump = std::max(m.d2_jump, std::abs(wrap_d2));
    }
    return m;
}

struct SecondDifferenceMetrics
{
    Real max_abs = 0.0;
    Real rms = 0.0;
};

SecondDifferenceMetrics secondDifferenceY(
    const amrex::Array4<const Real>& arr,
    const Box& region)
{
    SecondDifferenceMetrics m;
    Real sum_sq = 0.0;
    long count = 0;
    for (int j = region.smallEnd(1) + 1; j <= region.bigEnd(1) - 1; ++j) {
        for (int i = region.smallEnd(0); i <= region.bigEnd(0); ++i) {
            const Real d2 =
                arr(i, j + 1, 0) - 2.0 * arr(i, j, 0) + arr(i, j - 1, 0);
            m.max_abs = std::max(m.max_abs, std::abs(d2));
            sum_sq += d2 * d2;
            ++count;
        }
    }
    m.rms = (count > 0) ? std::sqrt(sum_sq / static_cast<Real>(count)) : 0.0;
    return m;
}

Real maxAbsSecondDifferenceX(
    const amrex::Array4<const Real>& arr,
    const Box& region)
{
    Real max_abs = 0.0;
    for (int j = region.smallEnd(1); j <= region.bigEnd(1); ++j) {
        for (int i = region.smallEnd(0) + 1; i <= region.bigEnd(0) - 1; ++i) {
            const Real d2 =
                arr(i + 1, j, 0) - 2.0 * arr(i, j, 0) + arr(i - 1, j, 0);
            max_abs = std::max(max_abs, std::abs(d2));
        }
    }
    return max_abs;
}

} // namespace

void writeNestedOscillationDiagnostics(
    crescent::shared::LevelContext& level_context,
    amrex::Real time_next)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    auto& runtime = level_context.runtime();
    auto& reporting = runtime.reporting;
    const auto& coupling_zones = level_context.nestedRuntime().coupling_zones;
    auto& periodic_fields =
        level_context.nestedRuntime().periodic_derivative_scratch[lev];
    auto& gather_root =
        level_context.nestedRuntime().runtime_diagnostic_scratch[lev].field_root[1];

    if (!crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev) ||
        reporting.write_omega_a_edge_diagnostics != 1) {
        return;
    }
    if (runtime.istep[lev] % reporting.omega_a_edge_diagnostic_int != 0) {
        return;
    }

    const auto& domain = coupling_zones.domains[lev];
    const Box& box_A = domain.box_A;

    // Sequential gathers reuse one rank-0 scratch: gather one field,
    // reduce its metrics on the IO rank, then gather the next.
    // The root scratch covers the whole level domain; the gathered data only
    // fills the L1 patch, so all index math must use the patch box.
    crescent::utils::gatherToRoot(solution.eta[lev], gather_root);
    const Box root_box = solution.eta[lev].boxArray().minimalBox();
    SecondDifferenceMetrics d2y_omega_a;
    Real d2x_omega_a_max = 0.0;
    SecondDifferenceMetrics d2y_band_lo;
    SecondDifferenceMetrics d2y_band_hi;
    if (ParallelDescriptor::IOProcessor()) {
        auto eta_arr = gather_root.const_arrays()[0];
        d2y_omega_a = secondDifferenceY(eta_arr, box_A);
        d2x_omega_a_max = maxAbsSecondDifferenceX(eta_arr, box_A);

        // 6-cell strips just outside the Omega_A y-edges (the documented
        // locus of the boundary packet), clamped to the level box.
        constexpr int band_width = 6;
        Box band_lo(
            amrex::IntVect(box_A.smallEnd(0), box_A.smallEnd(1) - band_width, 0),
            amrex::IntVect(box_A.bigEnd(0), box_A.smallEnd(1) - 1, 0));
        Box band_hi(
            amrex::IntVect(box_A.smallEnd(0), box_A.bigEnd(1) + 1, 0),
            amrex::IntVect(box_A.bigEnd(0), box_A.bigEnd(1) + band_width, 0));
        band_lo &= root_box;
        band_hi &= root_box;
        if (!band_lo.isEmpty()) {
            d2y_band_lo = secondDifferenceY(eta_arr, band_lo);
        }
        if (!band_hi.isEmpty()) {
            d2y_band_hi = secondDifferenceY(eta_arr, band_hi);
        }
    }

    // The seam columns monitor the fields that actually enter child FFTs:
    // the C3-periodized copies in the total-field modes, the tapered
    // correction pair in correction-field coupling (where the periodic
    // scratch is never written and must not be read).
    const bool correction_field_mode =
        level_context.nestedConfig().l1_coupling_mode == 2;
    const auto& correction_field =
        level_context.nestedRuntime().correction_field;
    const amrex::MultiFab& eta_seam_source = correction_field_mode
        ? correction_field.de
        : periodic_fields.eta_periodic;
    const amrex::MultiFab& phi_seam_source = correction_field_mode
        ? correction_field.dp
        : periodic_fields.phi_periodic;

    crescent::utils::gatherToRoot(eta_seam_source, gather_root);
    SeamMetrics eta_seam_x;
    SeamMetrics eta_seam_y;
    if (ParallelDescriptor::IOProcessor()) {
        auto eta_periodic_arr = gather_root.const_arrays()[0];
        eta_seam_x = seamMetricsX(eta_periodic_arr, root_box);
        eta_seam_y = seamMetricsY(eta_periodic_arr, root_box);
    }

    crescent::utils::gatherToRoot(phi_seam_source, gather_root);
    SeamMetrics phi_seam_x;
    SeamMetrics phi_seam_y;
    Real d2x_phi_omega_a_max = 0.0;
    SecondDifferenceMetrics d2y_phi_omega_a;
    if (ParallelDescriptor::IOProcessor()) {
        auto phi_periodic_arr = gather_root.const_arrays()[0];
        phi_seam_x = seamMetricsX(phi_periodic_arr, root_box);
        phi_seam_y = seamMetricsY(phi_periodic_arr, root_box);
        // Interior curvature of the same field: the reference scale that makes
        // the seam second-difference columns interpretable.
        d2x_phi_omega_a_max = maxAbsSecondDifferenceX(phi_periodic_arr, box_A);
        d2y_phi_omega_a = secondDifferenceY(phi_periodic_arr, box_A);
    }

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string out_dir = runtimeOutputDirectory(reporting.plot_file);
    const std::string filename = out_dir + "/nested_oscillation_diagnostic.csv";

    const bool first_write_this_run =
        !reporting.nested_oscillation_diagnostics_initialized ||
        !amrex::FileSystem::Exists(filename);
    const auto open_mode = first_write_this_run ? std::ios::trunc : std::ios::app;
    std::ofstream ofs(filename, open_mode);
    if (!ofs.is_open()) {
        amrex::Abort("Failed to open nested_oscillation_diagnostic.csv for writing");
    }
    if (first_write_this_run) {
        ofs << "step,time,"
               "max_d2y_eta_omegaA,rms_d2y_eta_omegaA,max_d2x_eta_omegaA,"
               "max_d2y_eta_yband_lo,max_d2y_eta_yband_hi,"
               "seamx_step_eta_periodic,seamx_d2_eta_periodic,"
               "seamy_step_eta_periodic,seamy_d2_eta_periodic,"
               "seamx_step_phi_periodic,seamx_d2_phi_periodic,"
               "seamy_step_phi_periodic,seamy_d2_phi_periodic,"
               "max_d2x_phi_omegaA,max_d2y_phi_omegaA\n";
        reporting.nested_oscillation_diagnostics_initialized = true;
    }

    ofs << runtime.istep[lev] << "," << time_next << ","
        << d2y_omega_a.max_abs << "," << d2y_omega_a.rms << ","
        << d2x_omega_a_max << ","
        << d2y_band_lo.max_abs << "," << d2y_band_hi.max_abs << ","
        << eta_seam_x.step_jump << "," << eta_seam_x.d2_jump << ","
        << eta_seam_y.step_jump << "," << eta_seam_y.d2_jump << ","
        << phi_seam_x.step_jump << "," << phi_seam_x.d2_jump << ","
        << phi_seam_y.step_jump << "," << phi_seam_y.d2_jump << ","
        << d2x_phi_omega_a_max << "," << d2y_phi_omega_a.max_abs << "\n";
}

} // namespace crescent::solver::runtime_diagnostics

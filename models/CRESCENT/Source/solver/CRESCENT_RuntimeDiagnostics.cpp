#include "solver/CRESCENT_RuntimeDiagnostics.H"

// Runtime Omega-A edge diagnostic implementation.

#include <fstream>
#include <vector>

#include <AMReX_FileSystem.H>

#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_Utils.H"
#include "solver/CRESCENT_RuntimeDiagnosticMetrics.H"

using namespace amrex;

namespace {
using namespace crescent::solver::runtime_diagnostics::detail;
}

namespace crescent::solver::runtime_diagnostics {

void writeOmegaAEdgeDiagnostics(
    crescent::shared::LevelEvolutionContext& level_context,
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
        level_context.nestedRuntime().runtime_diagnostic_scratch[lev].omega_a_edge;
    const auto& level_geom = level_context.levelGeometry();
    using namespace amrex;

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

    crescent::utils::gatherToRoot(solution.eta[lev], diagnostic_scratch.eta_root);
    crescent::utils::gatherToRoot(solution.phi[lev], diagnostic_scratch.phi_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.core_derivatives_on_l1.eta_x[lev],
        diagnostic_scratch.eta_x_core_interp_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.core_derivatives_on_l1.eta_y[lev],
        diagnostic_scratch.eta_y_core_interp_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.core_derivatives_on_l1.phi_x[lev],
        diagnostic_scratch.phi_x_core_interp_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.core_derivatives_on_l1.phi_y[lev],
        diagnostic_scratch.phi_y_core_interp_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.l1_blended_derivatives.eta_x[lev], diagnostic_scratch.eta_x_rk_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.l1_blended_derivatives.eta_y[lev], diagnostic_scratch.eta_y_rk_root);
    crescent::utils::gatherToRoot(solution.phi_x[lev], diagnostic_scratch.phi_x_root);
    crescent::utils::gatherToRoot(solution.phi_y[lev], diagnostic_scratch.phi_y_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.l1_blended_derivatives.phi_x[lev], diagnostic_scratch.phi_x_rk_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.l1_blended_derivatives.phi_y[lev], diagnostic_scratch.phi_y_rk_root);
    crescent::utils::gatherToRoot(
        coupling_derivatives.l1_blended_derivatives.phi_z[lev], diagnostic_scratch.phi_z_rk_root);
    crescent::utils::gatherToRoot(
        coupling_zones.sigma_blend[lev], diagnostic_scratch.sigma_root);

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string out_dir = runtimeOutputDirectory(reporting.plot_file);
    const std::string filename = out_dir + "/omega_a_edge_diagnostic.csv";
    const std::string spectral_filename = out_dir + "/omega_a_edge_spectral_diagnostic.csv";

    const bool first_write_this_run =
        !reporting.omega_a_edge_diagnostics_initialized ||
        !amrex::FileSystem::Exists(filename);
    const auto open_mode = first_write_this_run ? std::ios::trunc : std::ios::app;
    std::ofstream ofs(filename, open_mode);
    std::ofstream spectral_ofs(spectral_filename, open_mode);
    if (!ofs.is_open()) {
        amrex::Abort("Failed to open omega_a_edge_diagnostic.csv for writing");
    }
    if (!spectral_ofs.is_open()) {
        amrex::Abort("Failed to open omega_a_edge_spectral_diagnostic.csv for writing");
    }
    if (first_write_this_run) {
        ofs << "step,time,probe_label,i_probe,x_probe,j_lo,j_hi,"
               "sigma_min,sigma_max,eta_min,eta_max,phi_min,phi_max,"
               "eta_x_interp_min,eta_x_interp_max,eta_y_interp_min,eta_y_interp_max,"
               "phi_x_interp_min,phi_x_interp_max,phi_y_interp_min,phi_y_interp_max,"
               "eta_x_RK_min,eta_x_RK_max,eta_y_RK_min,eta_y_RK_max,"
               "phi_x_RK_min,phi_x_RK_max,phi_y_RK_min,phi_y_RK_max,"
               "phi_z_RK_min,phi_z_RK_max\n";
        spectral_ofs << "step,time,probe_label,i_probe,x_probe,field,n_samples,"
                        "mode0_energy,nonzero_energy,total_energy,nonzero_fraction,"
                        "mean_abs_mode,peak_abs_mode,peak_abs_mode_energy\n";
        reporting.omega_a_edge_diagnostics_initialized = true;
    }

    auto eta_arr = diagnostic_scratch.eta_root.const_arrays()[0];
    auto phi_arr = diagnostic_scratch.phi_root.const_arrays()[0];
    auto eta_x_core_interp_arr =
        diagnostic_scratch.eta_x_core_interp_root.const_arrays()[0];
    auto eta_y_core_interp_arr =
        diagnostic_scratch.eta_y_core_interp_root.const_arrays()[0];
    auto phi_x_core_interp_arr =
        diagnostic_scratch.phi_x_core_interp_root.const_arrays()[0];
    auto phi_y_core_interp_arr =
        diagnostic_scratch.phi_y_core_interp_root.const_arrays()[0];
    auto eta_x_RK_arr = diagnostic_scratch.eta_x_rk_root.const_arrays()[0];
    auto eta_y_RK_arr = diagnostic_scratch.eta_y_rk_root.const_arrays()[0];
    auto phi_x_arr = diagnostic_scratch.phi_x_root.const_arrays()[0];
    auto phi_y_arr = diagnostic_scratch.phi_y_root.const_arrays()[0];
    auto phi_x_RK_arr = diagnostic_scratch.phi_x_rk_root.const_arrays()[0];
    auto phi_y_RK_arr = diagnostic_scratch.phi_y_rk_root.const_arrays()[0];
    auto phi_z_RK_arr = diagnostic_scratch.phi_z_rk_root.const_arrays()[0];
    auto sigma_arr = diagnostic_scratch.sigma_root.const_arrays()[0];
    const Real* dx = level_geom.CellSize();
    const auto prob_lo = level_geom.ProbLoArray();

    auto write_line_spectrum = [&](const char* probe_label,
                                   int i_probe,
                                   Real x_probe,
                                   const char* field_name,
                                   const std::vector<Real>& samples)
    {
        const LineSpectrumMetrics metrics = summarizeLineSpectrum(samples);
        if (metrics.n <= 1) {
            return;
        }

        spectral_ofs << runtime.istep[lev] << "," << time_next << ","
                     << probe_label << "," << i_probe << "," << x_probe << ","
                     << field_name << "," << metrics.n << ","
                     << metrics.mode0_energy << "," << metrics.nonzero_energy << ","
                     << metrics.total_energy << "," << metrics.nonzero_fraction << ","
                     << metrics.mean_abs_mode << "," << metrics.peak_abs_mode << ","
                     << metrics.peak_abs_mode_energy << "\n";
    };

    for (int probe = 0; probe < 6; ++probe) {
        const int i_probe = i_candidates[probe];
        if (i_probe < level_domain.smallEnd(0) || i_probe > level_domain.bigEnd(0)) {
            continue;
        }

        Real sigma_min = 0.0, sigma_max = 0.0;
        Real eta_min = 0.0, eta_max = 0.0;
        Real phi_min = 0.0, phi_max = 0.0;
        Real eta_x_core_interp_min = 0.0, eta_x_core_interp_max = 0.0;
        Real eta_y_core_interp_min = 0.0, eta_y_core_interp_max = 0.0;
        Real phi_x_core_interp_min = 0.0, phi_x_core_interp_max = 0.0;
        Real phi_y_core_interp_min = 0.0, phi_y_core_interp_max = 0.0;
        Real eta_x_RK_min = 0.0, eta_x_RK_max = 0.0;
        Real eta_y_RK_min = 0.0, eta_y_RK_max = 0.0;
        Real phi_x_RK_min = 0.0, phi_x_RK_max = 0.0;
        Real phi_y_RK_min = 0.0, phi_y_RK_max = 0.0;
        Real phi_z_RK_min = 0.0, phi_z_RK_max = 0.0;
        std::vector<Real> phi_x_core_interp_line;
        std::vector<Real> phi_y_core_interp_line;
        std::vector<Real> phi_x_line;
        std::vector<Real> phi_y_line;
        std::vector<Real> phi_x_RK_line;
        std::vector<Real> phi_y_RK_line;
        phi_x_core_interp_line.reserve(j_hi - j_lo + 1);
        phi_y_core_interp_line.reserve(j_hi - j_lo + 1);
        phi_x_line.reserve(j_hi - j_lo + 1);
        phi_y_line.reserve(j_hi - j_lo + 1);
        phi_x_RK_line.reserve(j_hi - j_lo + 1);
        phi_y_RK_line.reserve(j_hi - j_lo + 1);

        for (int j = j_lo; j <= j_hi; ++j) {
            const Real sigma_val = sigma_arr(i_probe, j, 0);
            const Real eta_val = eta_arr(i_probe, j, 0);
            const Real phi_val = phi_arr(i_probe, j, 0);
            const Real eta_x_core_interp_val = eta_x_core_interp_arr(i_probe, j, 0);
            const Real eta_y_core_interp_val = eta_y_core_interp_arr(i_probe, j, 0);
            const Real phi_x_core_interp_val = phi_x_core_interp_arr(i_probe, j, 0);
            const Real phi_y_core_interp_val = phi_y_core_interp_arr(i_probe, j, 0);
            const Real eta_x_RK_val = eta_x_RK_arr(i_probe, j, 0);
            const Real eta_y_RK_val = eta_y_RK_arr(i_probe, j, 0);
            const Real phi_x_val = phi_x_arr(i_probe, j, 0);
            const Real phi_y_val = phi_y_arr(i_probe, j, 0);
            const Real phi_x_RK_val = phi_x_RK_arr(i_probe, j, 0);
            const Real phi_y_RK_val = phi_y_RK_arr(i_probe, j, 0);
            const Real phi_z_RK_val = phi_z_RK_arr(i_probe, j, 0);
            phi_x_core_interp_line.push_back(phi_x_core_interp_val);
            phi_y_core_interp_line.push_back(phi_y_core_interp_val);
            phi_x_line.push_back(phi_x_val);
            phi_y_line.push_back(phi_y_val);
            phi_x_RK_line.push_back(phi_x_RK_val);
            phi_y_RK_line.push_back(phi_y_RK_val);

            if (j == j_lo) {
                sigma_min = sigma_max = sigma_val;
                eta_min = eta_max = eta_val;
                phi_min = phi_max = phi_val;
                eta_x_core_interp_min = eta_x_core_interp_max = eta_x_core_interp_val;
                eta_y_core_interp_min = eta_y_core_interp_max = eta_y_core_interp_val;
                phi_x_core_interp_min = phi_x_core_interp_max = phi_x_core_interp_val;
                phi_y_core_interp_min = phi_y_core_interp_max = phi_y_core_interp_val;
                eta_x_RK_min = eta_x_RK_max = eta_x_RK_val;
                eta_y_RK_min = eta_y_RK_max = eta_y_RK_val;
                phi_x_RK_min = phi_x_RK_max = phi_x_RK_val;
                phi_y_RK_min = phi_y_RK_max = phi_y_RK_val;
                phi_z_RK_min = phi_z_RK_max = phi_z_RK_val;
            } else {
                sigma_min = std::min(sigma_min, sigma_val);
                sigma_max = std::max(sigma_max, sigma_val);
                eta_min = std::min(eta_min, eta_val);
                eta_max = std::max(eta_max, eta_val);
                phi_min = std::min(phi_min, phi_val);
                phi_max = std::max(phi_max, phi_val);
                eta_x_core_interp_min =
                    std::min(eta_x_core_interp_min, eta_x_core_interp_val);
                eta_x_core_interp_max =
                    std::max(eta_x_core_interp_max, eta_x_core_interp_val);
                eta_y_core_interp_min =
                    std::min(eta_y_core_interp_min, eta_y_core_interp_val);
                eta_y_core_interp_max =
                    std::max(eta_y_core_interp_max, eta_y_core_interp_val);
                phi_x_core_interp_min =
                    std::min(phi_x_core_interp_min, phi_x_core_interp_val);
                phi_x_core_interp_max =
                    std::max(phi_x_core_interp_max, phi_x_core_interp_val);
                phi_y_core_interp_min =
                    std::min(phi_y_core_interp_min, phi_y_core_interp_val);
                phi_y_core_interp_max =
                    std::max(phi_y_core_interp_max, phi_y_core_interp_val);
                eta_x_RK_min = std::min(eta_x_RK_min, eta_x_RK_val);
                eta_x_RK_max = std::max(eta_x_RK_max, eta_x_RK_val);
                eta_y_RK_min = std::min(eta_y_RK_min, eta_y_RK_val);
                eta_y_RK_max = std::max(eta_y_RK_max, eta_y_RK_val);
                phi_x_RK_min = std::min(phi_x_RK_min, phi_x_RK_val);
                phi_x_RK_max = std::max(phi_x_RK_max, phi_x_RK_val);
                phi_y_RK_min = std::min(phi_y_RK_min, phi_y_RK_val);
                phi_y_RK_max = std::max(phi_y_RK_max, phi_y_RK_val);
                phi_z_RK_min = std::min(phi_z_RK_min, phi_z_RK_val);
                phi_z_RK_max = std::max(phi_z_RK_max, phi_z_RK_val);
            }
        }

        const Real x_probe =
            prob_lo[0] + (static_cast<Real>(i_probe) + Real(0.5)) * dx[0];
        ofs << runtime.istep[lev] << "," << time_next << ","
            << labels[probe] << "," << i_probe << "," << x_probe << ","
            << j_lo << "," << j_hi << ","
            << sigma_min << "," << sigma_max << ","
            << eta_min << "," << eta_max << ","
            << phi_min << "," << phi_max << ","
            << eta_x_core_interp_min << "," << eta_x_core_interp_max << ","
            << eta_y_core_interp_min << "," << eta_y_core_interp_max << ","
            << phi_x_core_interp_min << "," << phi_x_core_interp_max << ","
            << phi_y_core_interp_min << "," << phi_y_core_interp_max << ","
            << eta_x_RK_min << "," << eta_x_RK_max << ","
            << eta_y_RK_min << "," << eta_y_RK_max << ","
            << phi_x_RK_min << "," << phi_x_RK_max << ","
            << phi_y_RK_min << "," << phi_y_RK_max << ","
            << phi_z_RK_min << "," << phi_z_RK_max << "\n";

        write_line_spectrum(
            labels[probe], i_probe, x_probe, "phi_x_interp", phi_x_core_interp_line);
        write_line_spectrum(
            labels[probe], i_probe, x_probe, "phi_y_interp", phi_y_core_interp_line);
        write_line_spectrum(labels[probe], i_probe, x_probe, "phi_x_L1", phi_x_line);
        write_line_spectrum(labels[probe], i_probe, x_probe, "phi_y_L1", phi_y_line);
        write_line_spectrum(labels[probe], i_probe, x_probe, "phi_x_RK", phi_x_RK_line);
        write_line_spectrum(labels[probe], i_probe, x_probe, "phi_y_RK", phi_y_RK_line);
    }
}

} // namespace crescent::solver::runtime_diagnostics

#include "solver/CRESCENT_RuntimeDiagnostics.H"

// Runtime stage-phi/interface diagnostic implementation.

#include <fstream>
#include <vector>

#include <AMReX_FileSystem.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedFreeSurfaceEvolution.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_Utils.H"
#include "CRESCENT_LevelEvolutionContext.H"
#include "solver/CRESCENT_RuntimeDiagnosticMetrics.H"

using namespace amrex;
using namespace crescent::constants;

namespace {
using namespace crescent::solver::runtime_diagnostics::detail;
}

namespace crescent::solver::runtime_diagnostics {

void writeOmegaAEdgePhiStageDiagnostics(
    const crescent::shared::LevelEvolutionContext& level_context,
    const int irk,
    amrex::Real time_stage,
    const amrex::MultiFab& phi_raw,
    const amrex::MultiFab& phi_c3,
    const amrex::MultiFab& phi_filtered)
{
    using namespace amrex;
    const int lev = level_context.level();
    auto& runtime = level_context.runtime();
    auto& reporting = runtime.reporting;
    auto& spectral_workspace = level_context.spectralWorkspace();
    const auto& nested_state = level_context.nestedRuntime();
    auto& diagnostic_scratch =
        level_context.nestedRuntime().runtime_diagnostic_scratch[lev].phi_stage;
    const auto& coupling_zones = nested_state.coupling_zones;
    const auto& level_geom = level_context.levelGeometry();

    if (!crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev) ||
        irk != 0 ||
        reporting.write_omega_a_edge_diagnostics != 1) {
        return;
    }
    if (runtime.istep[lev] % reporting.omega_a_edge_diagnostic_int != 0) {
        return;
    }

    const auto& domain = coupling_zones.domains[lev];
    const Box& box_A = domain.box_A;
    const Box& box_S = domain.box_S;
    const Box level_domain = phi_raw.boxArray().minimalBox();
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

    auto compute_phi_y_fft = [&](const MultiFab& phi_in, MultiFab& phi_y_out)
    {
        const Box& box_nested = phi_in.boxArray().minimalBox();
        const int ni_nested = box_nested.length(0);
        const int nj_nested = box_nested.length(1);
        const Real dx_local = level_geom.CellSize(0);
        const Real dy_local = level_geom.CellSize(1);
        const Real dk_y_nested = PI2 / (static_cast<Real>(nj_nested) * dy_local);
        const int cutoff_x = static_cast<int>(ni_nested * ONE_THIRD_CUTOFF);
        const int cutoff_y = static_cast<int>(nj_nested * ONE_THIRD_CUTOFF);

        // Reuse the per-level directional derivative spectra scratch. This
        // diagnostic runs after stage derivative preparation returns, so the
        // buffers are free for sequential reuse here.
        auto& phi_c = spectral_workspace.directional_input_spectral[lev];
        auto& phi_y_c = spectral_workspace.directional_y_spectral[lev];

        crescent::spectral::forwardFFT(
            spectral_workspace.fft[lev], level_geom.data(), phi_in, phi_c);

        for (MFIter mfi(phi_c, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto phi_c_arr = phi_c.const_array(mfi);
            auto phi_y_c_arr = phi_y_c.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const int kx = (i <= ni_nested/2) ? i : i - ni_nested;
                const int ky_index = (j <= nj_nested/2) ? j : (j - nj_nested);
                const bool dealias_cutoff =
                    (std::abs(kx) > cutoff_x || std::abs(ky_index) > cutoff_y);

                if (dealias_cutoff) {
                    phi_y_c_arr(i,j,k) = GpuComplex<Real>(0.0, 0.0);
                } else {
                    const Real norm = Real(1.0) /
                                      (static_cast<Real>(ni_nested) *
                                       static_cast<Real>(nj_nested));
                    GpuComplex<Real> phase_y(0.0, ky_index * dk_y_nested);
                    phi_y_c_arr(i,j,k) = phi_c_arr(i,j,k) * phase_y * norm;
                }
            });
        }

        crescent::spectral::backwardFFT(
            spectral_workspace.fft[lev], level_geom.data(), phi_y_c, phi_y_out);
    };

    compute_phi_y_fft(phi_raw, diagnostic_scratch.phi_y_raw);
    compute_phi_y_fft(phi_c3, diagnostic_scratch.phi_y_c3);
    compute_phi_y_fft(phi_filtered, diagnostic_scratch.phi_y_filtered);

    crescent::utils::gatherToRoot(phi_raw, diagnostic_scratch.phi_raw_root);
    crescent::utils::gatherToRoot(phi_c3, diagnostic_scratch.phi_c3_root);
    crescent::utils::gatherToRoot(phi_filtered, diagnostic_scratch.phi_filtered_root);
    crescent::utils::gatherToRoot(
        diagnostic_scratch.phi_y_raw, diagnostic_scratch.phi_y_raw_root);
    crescent::utils::gatherToRoot(
        diagnostic_scratch.phi_y_c3, diagnostic_scratch.phi_y_c3_root);
    crescent::utils::gatherToRoot(
        diagnostic_scratch.phi_y_filtered, diagnostic_scratch.phi_y_filtered_root);

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string out_dir = runtimeOutputDirectory(reporting.plot_file);
    const std::string filename = out_dir + "/omega_a_edge_phi_stage_spectral_diagnostic.csv";
    const std::string interface_spectral_filename = out_dir + "/omega_s_y_interface_spectral_diagnostic.csv";
    const std::string interface_jump_filename = out_dir + "/omega_s_y_interface_jump_diagnostic.csv";
    const std::string x_interface_spectral_filename = out_dir + "/omega_s_x_interface_spectral_diagnostic.csv";
    const std::string x_interface_jump_filename = out_dir + "/omega_s_x_interface_jump_diagnostic.csv";
    const std::string corner_filename = out_dir + "/omega_s_corner_diagnostic.csv";

    const bool first_write_this_run =
        !reporting.omega_a_stage_diagnostics_initialized ||
        !amrex::FileSystem::Exists(filename);
    const auto open_mode = first_write_this_run ? std::ios::trunc : std::ios::app;
    std::ofstream ofs(filename, open_mode);
    std::ofstream interface_spectral_ofs(interface_spectral_filename, open_mode);
    std::ofstream interface_jump_ofs(interface_jump_filename, open_mode);
    std::ofstream x_interface_spectral_ofs(x_interface_spectral_filename, open_mode);
    std::ofstream x_interface_jump_ofs(x_interface_jump_filename, open_mode);
    std::ofstream corner_ofs(corner_filename, open_mode);
    if (!ofs.is_open()) {
        amrex::Abort("Failed to open omega_a_edge_phi_stage_spectral_diagnostic.csv for writing");
    }
    if (!interface_spectral_ofs.is_open()) {
        amrex::Abort("Failed to open omega_s_y_interface_spectral_diagnostic.csv for writing");
    }
    if (!interface_jump_ofs.is_open()) {
        amrex::Abort("Failed to open omega_s_y_interface_jump_diagnostic.csv for writing");
    }
    if (!x_interface_spectral_ofs.is_open()) {
        amrex::Abort("Failed to open omega_s_x_interface_spectral_diagnostic.csv for writing");
    }
    if (!x_interface_jump_ofs.is_open()) {
        amrex::Abort("Failed to open omega_s_x_interface_jump_diagnostic.csv for writing");
    }
    if (!corner_ofs.is_open()) {
        amrex::Abort("Failed to open omega_s_corner_diagnostic.csv for writing");
    }
    if (first_write_this_run) {
        ofs << "step,time,irk,probe_label,i_probe,x_probe,field,n_samples,"
               "mode0_energy,nonzero_energy,total_energy,nonzero_fraction,"
               "mean_abs_mode,peak_abs_mode,peak_abs_mode_energy\n";
        interface_spectral_ofs << "step,time,irk,row_label,j_probe,y_probe,field,n_samples,"
                                  "mode0_energy,nonzero_energy,total_energy,nonzero_fraction,"
                                  "mean_abs_mode,peak_abs_mode,peak_abs_mode_energy\n";
        interface_jump_ofs << "step,time,irk,interface_label,field,j0,j1,y0,y1,n_samples,"
                              "min_jump,max_jump,max_abs_jump,mean_abs_jump,rms_jump\n";
        x_interface_spectral_ofs << "step,time,irk,col_label,i_probe,x_probe,field,n_samples,"
                                    "mode0_energy,nonzero_energy,total_energy,nonzero_fraction,"
                                    "mean_abs_mode,peak_abs_mode,peak_abs_mode_energy\n";
        x_interface_jump_ofs << "step,time,irk,interface_label,field,i0,i1,x0,x1,n_samples,"
                                "min_jump,max_jump,max_abs_jump,mean_abs_jump,rms_jump\n";
        corner_ofs << "step,time,irk,corner_label,field,i_lo,i_hi,j_lo,j_hi,nx,ny,"
                      "min_val,max_val,max_abs_val,range_val,rms_val,max_abs_dx,max_abs_dy\n";
        reporting.omega_a_stage_diagnostics_initialized = true;
    }

    auto phi_raw_arr = diagnostic_scratch.phi_raw_root.const_arrays()[0];
    auto phi_c3_arr = diagnostic_scratch.phi_c3_root.const_arrays()[0];
    auto phi_filtered_arr = diagnostic_scratch.phi_filtered_root.const_arrays()[0];
    auto phi_y_raw_arr = diagnostic_scratch.phi_y_raw_root.const_arrays()[0];
    auto phi_y_c3_arr = diagnostic_scratch.phi_y_c3_root.const_arrays()[0];
    auto phi_y_filtered_arr =
        diagnostic_scratch.phi_y_filtered_root.const_arrays()[0];
    const Real* dx = level_geom.CellSize();
    const auto prob_lo = level_geom.ProbLoArray();

    auto write_line_spectrum = [&](const char* probe_label,
                                   int i_probe,
                                   Real x_probe,
                                   const char* field_name,
                                   const std::vector<Real>& samples,
                                   std::ofstream& stream)
    {
        const LineSpectrumMetrics metrics = summarizeLineSpectrum(samples);
        if (metrics.n <= 1) {
            return;
        }

        stream << runtime.istep[lev] << "," << time_stage << "," << irk << ","
               << probe_label << "," << i_probe << "," << x_probe << ","
               << field_name << "," << metrics.n << ","
               << metrics.mode0_energy << "," << metrics.nonzero_energy << ","
               << metrics.total_energy << "," << metrics.nonzero_fraction << ","
               << metrics.mean_abs_mode << "," << metrics.peak_abs_mode << ","
               << metrics.peak_abs_mode_energy << "\n";
    };

    auto write_row_jump = [&](const char* interface_label,
                              const char* field_name,
                              int j0,
                              int j1,
                              Real y0,
                              Real y1,
                              const std::vector<Real>& row0,
                              const std::vector<Real>& row1)
    {
        const JumpMetrics metrics = summarizeJumpMetrics(row0, row1);
        if (metrics.n <= 0) {
            return;
        }

        interface_jump_ofs << runtime.istep[lev] << "," << time_stage << "," << irk << ","
                           << interface_label << "," << field_name << ","
                           << j0 << "," << j1 << "," << y0 << "," << y1 << ","
                           << metrics.n << ","
                           << metrics.min_jump << "," << metrics.max_jump << ","
                           << metrics.max_abs_jump << "," << metrics.mean_abs_jump
                           << "," << metrics.rms_jump << "\n";
    };

    auto write_col_jump = [&](const char* interface_label,
                              const char* field_name,
                              int i0,
                              int i1,
                              Real x0,
                              Real x1,
                              const std::vector<Real>& col0,
                              const std::vector<Real>& col1)
    {
        const JumpMetrics metrics = summarizeJumpMetrics(col0, col1);
        if (metrics.n <= 0) {
            return;
        }

        x_interface_jump_ofs << runtime.istep[lev] << "," << time_stage << "," << irk << ","
                             << interface_label << "," << field_name << ","
                             << i0 << "," << i1 << "," << x0 << "," << x1 << ","
                             << metrics.n << ","
                             << metrics.min_jump << "," << metrics.max_jump << ","
                             << metrics.max_abs_jump << "," << metrics.mean_abs_jump
                             << "," << metrics.rms_jump << "\n";
    };

    auto write_patch_stats = [&](const char* corner_label,
                                 const char* field_name,
                                 int i_lo_patch,
                                 int i_hi_patch,
                                 int j_lo_patch,
                                 int j_hi_patch,
                                 auto const& arr)
    {
        PatchMetrics metrics;
        if (!computePatchMetrics(i_lo_patch, i_hi_patch, j_lo_patch, j_hi_patch, arr, metrics)) {
            return;
        }

        corner_ofs << runtime.istep[lev] << "," << time_stage << "," << irk << ","
                   << corner_label << "," << field_name << ","
                   << i_lo_patch << "," << i_hi_patch << ","
                   << j_lo_patch << "," << j_hi_patch << ","
                   << metrics.nx_patch << "," << metrics.ny_patch << ","
                   << metrics.min_val << "," << metrics.max_val << ","
                   << metrics.max_abs_val << ","
                   << (metrics.max_val - metrics.min_val) << ","
                   << metrics.rms_val << ","
                   << metrics.max_abs_dx << "," << metrics.max_abs_dy << "\n";
    };

    for (int probe = 0; probe < 6; ++probe) {
        const int i_probe = i_candidates[probe];
        if (i_probe < level_domain.smallEnd(0) || i_probe > level_domain.bigEnd(0)) {
            continue;
        }

        std::vector<Real> phi_raw_line;
        std::vector<Real> phi_c3_line;
        std::vector<Real> phi_filtered_line;
        std::vector<Real> phi_y_raw_line;
        std::vector<Real> phi_y_c3_line;
        std::vector<Real> phi_y_filtered_line;
        phi_raw_line.reserve(j_hi - j_lo + 1);
        phi_c3_line.reserve(j_hi - j_lo + 1);
        phi_filtered_line.reserve(j_hi - j_lo + 1);
        phi_y_raw_line.reserve(j_hi - j_lo + 1);
        phi_y_c3_line.reserve(j_hi - j_lo + 1);
        phi_y_filtered_line.reserve(j_hi - j_lo + 1);

        for (int j = j_lo; j <= j_hi; ++j) {
            phi_raw_line.push_back(phi_raw_arr(i_probe, j, 0));
            phi_c3_line.push_back(phi_c3_arr(i_probe, j, 0));
            phi_filtered_line.push_back(phi_filtered_arr(i_probe, j, 0));
            phi_y_raw_line.push_back(phi_y_raw_arr(i_probe, j, 0));
            phi_y_c3_line.push_back(phi_y_c3_arr(i_probe, j, 0));
            phi_y_filtered_line.push_back(phi_y_filtered_arr(i_probe, j, 0));
        }

        const Real x_probe =
            prob_lo[0] + (static_cast<Real>(i_probe) + Real(0.5)) * dx[0];
        write_line_spectrum(labels[probe], i_probe, x_probe, "phi_raw", phi_raw_line, ofs);
        write_line_spectrum(labels[probe], i_probe, x_probe, "phi_c3", phi_c3_line, ofs);
        write_line_spectrum(
            labels[probe], i_probe, x_probe, "phi_filtered", phi_filtered_line, ofs);
        write_line_spectrum(
            labels[probe], i_probe, x_probe, "phi_y_rawfft", phi_y_raw_line, ofs);
        write_line_spectrum(
            labels[probe], i_probe, x_probe, "phi_y_c3fft", phi_y_c3_line, ofs);
        write_line_spectrum(
            labels[probe], i_probe, x_probe, "phi_y_filteredfft", phi_y_filtered_line, ofs);
    }

    const int i_lo = level_domain.smallEnd(0);
    const int i_hi = level_domain.bigEnd(0);
    const int j_bottom = box_S.smallEnd(1);
    const int j_top = box_S.bigEnd(1);
    const int row_candidates[] = {
        j_bottom - 1,
        j_bottom,
        j_bottom + 1,
        j_top - 1,
        j_top,
        j_top + 1
    };
    const char* row_labels[] = {
        "bottom_outer",
        "bottom_interface",
        "bottom_inner",
        "top_inner",
        "top_interface",
        "top_outer"
    };

    const int j_min = level_domain.smallEnd(1);
    const int j_max = level_domain.bigEnd(1);

    for (int row = 0; row < 6; ++row) {
        const int j_probe = row_candidates[row];
        if (j_probe < level_domain.smallEnd(1) || j_probe > level_domain.bigEnd(1)) {
            continue;
        }

        const Real y_probe =
            prob_lo[1] + (static_cast<Real>(j_probe) + Real(0.5)) * dx[1];
        const auto phi_raw_row = extractRowSamples(phi_raw_arr, i_lo, i_hi, j_probe);
        const auto phi_c3_row = extractRowSamples(phi_c3_arr, i_lo, i_hi, j_probe);
        const auto phi_filtered_row =
            extractRowSamples(phi_filtered_arr, i_lo, i_hi, j_probe);
        const auto phi_y_raw_row =
            extractRowSamples(phi_y_raw_arr, i_lo, i_hi, j_probe);
        const auto phi_y_c3_row =
            extractRowSamples(phi_y_c3_arr, i_lo, i_hi, j_probe);
        const auto phi_y_filtered_row =
            extractRowSamples(phi_y_filtered_arr, i_lo, i_hi, j_probe);

        write_line_spectrum(
            row_labels[row], j_probe, y_probe, "phi_raw", phi_raw_row, interface_spectral_ofs);
        write_line_spectrum(
            row_labels[row], j_probe, y_probe, "phi_c3", phi_c3_row, interface_spectral_ofs);
        write_line_spectrum(
            row_labels[row], j_probe, y_probe, "phi_filtered", phi_filtered_row,
            interface_spectral_ofs);
        write_line_spectrum(
            row_labels[row], j_probe, y_probe, "phi_y_rawfft", phi_y_raw_row,
            interface_spectral_ofs);
        write_line_spectrum(
            row_labels[row], j_probe, y_probe, "phi_y_c3fft", phi_y_c3_row,
            interface_spectral_ofs);
        write_line_spectrum(
            row_labels[row], j_probe, y_probe, "phi_y_filteredfft", phi_y_filtered_row,
            interface_spectral_ofs);
    }

    if (j_bottom - 1 >= level_domain.smallEnd(1) && j_bottom <= level_domain.bigEnd(1)) {
        const Real y0 =
            prob_lo[1] + (static_cast<Real>(j_bottom - 1) + Real(0.5)) * dx[1];
        const Real y1 =
            prob_lo[1] + (static_cast<Real>(j_bottom) + Real(0.5)) * dx[1];
        const auto phi_raw_outer =
            extractRowSamples(phi_raw_arr, i_lo, i_hi, j_bottom - 1);
        const auto phi_raw_inner = extractRowSamples(phi_raw_arr, i_lo, i_hi, j_bottom);
        const auto phi_c3_outer =
            extractRowSamples(phi_c3_arr, i_lo, i_hi, j_bottom - 1);
        const auto phi_c3_inner = extractRowSamples(phi_c3_arr, i_lo, i_hi, j_bottom);
        const auto phi_filtered_outer =
            extractRowSamples(phi_filtered_arr, i_lo, i_hi, j_bottom - 1);
        const auto phi_filtered_inner =
            extractRowSamples(phi_filtered_arr, i_lo, i_hi, j_bottom);
        const auto phi_y_raw_outer =
            extractRowSamples(phi_y_raw_arr, i_lo, i_hi, j_bottom - 1);
        const auto phi_y_raw_inner =
            extractRowSamples(phi_y_raw_arr, i_lo, i_hi, j_bottom);
        const auto phi_y_c3_outer =
            extractRowSamples(phi_y_c3_arr, i_lo, i_hi, j_bottom - 1);
        const auto phi_y_c3_inner =
            extractRowSamples(phi_y_c3_arr, i_lo, i_hi, j_bottom);
        const auto phi_y_filtered_outer =
            extractRowSamples(phi_y_filtered_arr, i_lo, i_hi, j_bottom - 1);
        const auto phi_y_filtered_inner =
            extractRowSamples(phi_y_filtered_arr, i_lo, i_hi, j_bottom);

        write_row_jump(
            "bottom_interface", "phi_raw", j_bottom - 1, j_bottom, y0, y1,
            phi_raw_outer, phi_raw_inner);
        write_row_jump(
            "bottom_interface", "phi_c3", j_bottom - 1, j_bottom, y0, y1,
            phi_c3_outer, phi_c3_inner);
        write_row_jump(
            "bottom_interface", "phi_filtered", j_bottom - 1, j_bottom, y0, y1,
            phi_filtered_outer, phi_filtered_inner);
        write_row_jump(
            "bottom_interface", "phi_y_rawfft", j_bottom - 1, j_bottom, y0, y1,
            phi_y_raw_outer, phi_y_raw_inner);
        write_row_jump(
            "bottom_interface", "phi_y_c3fft", j_bottom - 1, j_bottom, y0, y1,
            phi_y_c3_outer, phi_y_c3_inner);
        write_row_jump(
            "bottom_interface", "phi_y_filteredfft", j_bottom - 1, j_bottom, y0, y1,
            phi_y_filtered_outer, phi_y_filtered_inner);
    }

    if (j_top >= level_domain.smallEnd(1) && j_top + 1 <= level_domain.bigEnd(1)) {
        const Real y0 =
            prob_lo[1] + (static_cast<Real>(j_top) + Real(0.5)) * dx[1];
        const Real y1 =
            prob_lo[1] + (static_cast<Real>(j_top + 1) + Real(0.5)) * dx[1];
        const auto phi_raw_inner = extractRowSamples(phi_raw_arr, i_lo, i_hi, j_top);
        const auto phi_raw_outer =
            extractRowSamples(phi_raw_arr, i_lo, i_hi, j_top + 1);
        const auto phi_c3_inner = extractRowSamples(phi_c3_arr, i_lo, i_hi, j_top);
        const auto phi_c3_outer =
            extractRowSamples(phi_c3_arr, i_lo, i_hi, j_top + 1);
        const auto phi_filtered_inner =
            extractRowSamples(phi_filtered_arr, i_lo, i_hi, j_top);
        const auto phi_filtered_outer =
            extractRowSamples(phi_filtered_arr, i_lo, i_hi, j_top + 1);
        const auto phi_y_raw_inner =
            extractRowSamples(phi_y_raw_arr, i_lo, i_hi, j_top);
        const auto phi_y_raw_outer =
            extractRowSamples(phi_y_raw_arr, i_lo, i_hi, j_top + 1);
        const auto phi_y_c3_inner =
            extractRowSamples(phi_y_c3_arr, i_lo, i_hi, j_top);
        const auto phi_y_c3_outer =
            extractRowSamples(phi_y_c3_arr, i_lo, i_hi, j_top + 1);
        const auto phi_y_filtered_inner =
            extractRowSamples(phi_y_filtered_arr, i_lo, i_hi, j_top);
        const auto phi_y_filtered_outer =
            extractRowSamples(phi_y_filtered_arr, i_lo, i_hi, j_top + 1);

        write_row_jump(
            "top_interface", "phi_raw", j_top, j_top + 1, y0, y1,
            phi_raw_inner, phi_raw_outer);
        write_row_jump(
            "top_interface", "phi_c3", j_top, j_top + 1, y0, y1,
            phi_c3_inner, phi_c3_outer);
        write_row_jump(
            "top_interface", "phi_filtered", j_top, j_top + 1, y0, y1,
            phi_filtered_inner, phi_filtered_outer);
        write_row_jump(
            "top_interface", "phi_y_rawfft", j_top, j_top + 1, y0, y1,
            phi_y_raw_inner, phi_y_raw_outer);
        write_row_jump(
            "top_interface", "phi_y_c3fft", j_top, j_top + 1, y0, y1,
            phi_y_c3_inner, phi_y_c3_outer);
        write_row_jump(
            "top_interface", "phi_y_filteredfft", j_top, j_top + 1, y0, y1,
            phi_y_filtered_inner, phi_y_filtered_outer);
    }

    const int i_left = box_S.smallEnd(0);
    const int i_right = box_S.bigEnd(0);
    const int col_candidates[] = {
        i_left - 1,
        i_left,
        i_left + 1,
        i_right - 1,
        i_right,
        i_right + 1
    };
    const char* col_labels[] = {
        "left_outer",
        "left_interface",
        "left_inner",
        "right_inner",
        "right_interface",
        "right_outer"
    };

    for (int col = 0; col < 6; ++col) {
        const int i_probe = col_candidates[col];
        if (i_probe < level_domain.smallEnd(0) || i_probe > level_domain.bigEnd(0)) {
            continue;
        }

        const Real x_probe =
            prob_lo[0] + (static_cast<Real>(i_probe) + Real(0.5)) * dx[0];
        const auto phi_raw_col = extractColumnSamples(phi_raw_arr, j_min, j_max, i_probe);
        const auto phi_c3_col = extractColumnSamples(phi_c3_arr, j_min, j_max, i_probe);
        const auto phi_filtered_col =
            extractColumnSamples(phi_filtered_arr, j_min, j_max, i_probe);
        const auto phi_y_raw_col =
            extractColumnSamples(phi_y_raw_arr, j_min, j_max, i_probe);
        const auto phi_y_c3_col =
            extractColumnSamples(phi_y_c3_arr, j_min, j_max, i_probe);
        const auto phi_y_filtered_col =
            extractColumnSamples(phi_y_filtered_arr, j_min, j_max, i_probe);

        write_line_spectrum(
            col_labels[col], i_probe, x_probe, "phi_raw", phi_raw_col, x_interface_spectral_ofs);
        write_line_spectrum(
            col_labels[col], i_probe, x_probe, "phi_c3", phi_c3_col, x_interface_spectral_ofs);
        write_line_spectrum(
            col_labels[col], i_probe, x_probe, "phi_filtered", phi_filtered_col,
            x_interface_spectral_ofs);
        write_line_spectrum(
            col_labels[col], i_probe, x_probe, "phi_y_rawfft", phi_y_raw_col,
            x_interface_spectral_ofs);
        write_line_spectrum(
            col_labels[col], i_probe, x_probe, "phi_y_c3fft", phi_y_c3_col,
            x_interface_spectral_ofs);
        write_line_spectrum(
            col_labels[col], i_probe, x_probe, "phi_y_filteredfft", phi_y_filtered_col,
            x_interface_spectral_ofs);
    }

    if (i_left - 1 >= level_domain.smallEnd(0) && i_left <= level_domain.bigEnd(0)) {
        const Real x0 =
            prob_lo[0] + (static_cast<Real>(i_left - 1) + Real(0.5)) * dx[0];
        const Real x1 =
            prob_lo[0] + (static_cast<Real>(i_left) + Real(0.5)) * dx[0];
        const auto phi_raw_outer =
            extractColumnSamples(phi_raw_arr, j_min, j_max, i_left - 1);
        const auto phi_raw_inner =
            extractColumnSamples(phi_raw_arr, j_min, j_max, i_left);
        const auto phi_c3_outer =
            extractColumnSamples(phi_c3_arr, j_min, j_max, i_left - 1);
        const auto phi_c3_inner =
            extractColumnSamples(phi_c3_arr, j_min, j_max, i_left);
        const auto phi_filtered_outer =
            extractColumnSamples(phi_filtered_arr, j_min, j_max, i_left - 1);
        const auto phi_filtered_inner =
            extractColumnSamples(phi_filtered_arr, j_min, j_max, i_left);
        const auto phi_y_raw_outer =
            extractColumnSamples(phi_y_raw_arr, j_min, j_max, i_left - 1);
        const auto phi_y_raw_inner =
            extractColumnSamples(phi_y_raw_arr, j_min, j_max, i_left);
        const auto phi_y_c3_outer =
            extractColumnSamples(phi_y_c3_arr, j_min, j_max, i_left - 1);
        const auto phi_y_c3_inner =
            extractColumnSamples(phi_y_c3_arr, j_min, j_max, i_left);
        const auto phi_y_filtered_outer =
            extractColumnSamples(phi_y_filtered_arr, j_min, j_max, i_left - 1);
        const auto phi_y_filtered_inner =
            extractColumnSamples(phi_y_filtered_arr, j_min, j_max, i_left);

        write_col_jump(
            "left_interface", "phi_raw", i_left - 1, i_left, x0, x1,
            phi_raw_outer, phi_raw_inner);
        write_col_jump(
            "left_interface", "phi_c3", i_left - 1, i_left, x0, x1,
            phi_c3_outer, phi_c3_inner);
        write_col_jump(
            "left_interface", "phi_filtered", i_left - 1, i_left, x0, x1,
            phi_filtered_outer, phi_filtered_inner);
        write_col_jump(
            "left_interface", "phi_y_rawfft", i_left - 1, i_left, x0, x1,
            phi_y_raw_outer, phi_y_raw_inner);
        write_col_jump(
            "left_interface", "phi_y_c3fft", i_left - 1, i_left, x0, x1,
            phi_y_c3_outer, phi_y_c3_inner);
        write_col_jump(
            "left_interface", "phi_y_filteredfft", i_left - 1, i_left, x0, x1,
            phi_y_filtered_outer, phi_y_filtered_inner);
    }

    if (i_right >= level_domain.smallEnd(0) && i_right + 1 <= level_domain.bigEnd(0)) {
        const Real x0 =
            prob_lo[0] + (static_cast<Real>(i_right) + Real(0.5)) * dx[0];
        const Real x1 =
            prob_lo[0] + (static_cast<Real>(i_right + 1) + Real(0.5)) * dx[0];
        const auto phi_raw_inner =
            extractColumnSamples(phi_raw_arr, j_min, j_max, i_right);
        const auto phi_raw_outer =
            extractColumnSamples(phi_raw_arr, j_min, j_max, i_right + 1);
        const auto phi_c3_inner =
            extractColumnSamples(phi_c3_arr, j_min, j_max, i_right);
        const auto phi_c3_outer =
            extractColumnSamples(phi_c3_arr, j_min, j_max, i_right + 1);
        const auto phi_filtered_inner =
            extractColumnSamples(phi_filtered_arr, j_min, j_max, i_right);
        const auto phi_filtered_outer =
            extractColumnSamples(phi_filtered_arr, j_min, j_max, i_right + 1);
        const auto phi_y_raw_inner =
            extractColumnSamples(phi_y_raw_arr, j_min, j_max, i_right);
        const auto phi_y_raw_outer =
            extractColumnSamples(phi_y_raw_arr, j_min, j_max, i_right + 1);
        const auto phi_y_c3_inner =
            extractColumnSamples(phi_y_c3_arr, j_min, j_max, i_right);
        const auto phi_y_c3_outer =
            extractColumnSamples(phi_y_c3_arr, j_min, j_max, i_right + 1);
        const auto phi_y_filtered_inner =
            extractColumnSamples(phi_y_filtered_arr, j_min, j_max, i_right);
        const auto phi_y_filtered_outer =
            extractColumnSamples(phi_y_filtered_arr, j_min, j_max, i_right + 1);

        write_col_jump(
            "right_interface", "phi_raw", i_right, i_right + 1, x0, x1,
            phi_raw_inner, phi_raw_outer);
        write_col_jump(
            "right_interface", "phi_c3", i_right, i_right + 1, x0, x1,
            phi_c3_inner, phi_c3_outer);
        write_col_jump(
            "right_interface", "phi_filtered", i_right, i_right + 1, x0, x1,
            phi_filtered_inner, phi_filtered_outer);
        write_col_jump(
            "right_interface", "phi_y_rawfft", i_right, i_right + 1, x0, x1,
            phi_y_raw_inner, phi_y_raw_outer);
        write_col_jump(
            "right_interface", "phi_y_c3fft", i_right, i_right + 1, x0, x1,
            phi_y_c3_inner, phi_y_c3_outer);
        write_col_jump(
            "right_interface", "phi_y_filteredfft", i_right, i_right + 1, x0, x1,
            phi_y_filtered_inner, phi_y_filtered_outer);
    }

    const int corner_radius = 2;
    const struct CornerSpec {
        const char* label;
        int i_center;
        int j_center;
    } corners[] = {
        {"bottom_left",  i_left,  j_bottom},
        {"bottom_right", i_right, j_bottom},
        {"top_left",     i_left,  j_top},
        {"top_right",    i_right, j_top}
    };

    for (const auto& corner : corners) {
        const int i_lo_patch =
            std::max(level_domain.smallEnd(0), corner.i_center - corner_radius);
        const int i_hi_patch =
            std::min(level_domain.bigEnd(0), corner.i_center + corner_radius);
        const int j_lo_patch =
            std::max(level_domain.smallEnd(1), corner.j_center - corner_radius);
        const int j_hi_patch =
            std::min(level_domain.bigEnd(1), corner.j_center + corner_radius);

        write_patch_stats(
            corner.label, "phi_c3", i_lo_patch, i_hi_patch, j_lo_patch, j_hi_patch,
            phi_c3_arr);
        write_patch_stats(
            corner.label, "phi_filtered", i_lo_patch, i_hi_patch, j_lo_patch, j_hi_patch,
            phi_filtered_arr);
        write_patch_stats(
            corner.label, "phi_y_c3fft", i_lo_patch, i_hi_patch, j_lo_patch, j_hi_patch,
            phi_y_c3_arr);
        write_patch_stats(
            corner.label, "phi_y_filteredfft", i_lo_patch, i_hi_patch, j_lo_patch,
            j_hi_patch, phi_y_filtered_arr);
    }
}

} // namespace crescent::solver::runtime_diagnostics

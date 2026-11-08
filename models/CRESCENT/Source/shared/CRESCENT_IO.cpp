#include "CRESCENT_IO.H"
#include "CRESCENT_NestedDomains.H"

#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace amrex;

namespace crescent::io {

namespace {
    void openOutputFileOrAbort(
        std::ofstream& output,
        const std::string& filename,
        const char* description)
    {
        output.open(filename);
        if (!output.is_open()) {
            amrex::Abort(
                std::string("Failed to open ") + description + " '" +
                filename + "' for writing");
        }
    }
}

void writePlotfile(const std::string& plot_file,
                   const int step,
                   const int finest_level,
                   const Vector<std::string>& varnames,
                   const Vector<MultiFab>& mf,
                   const Vector<Geometry>& geom,
                   const Real time,
                   const Vector<int>& istep,
                   const Vector<IntVect>& ref_ratio)
{
    std::string plotfilename = Concatenate(plot_file, step, 5);

    const int mf_out_ncomp = varnames.size();
    if (mf_out_ncomp == 0) return;

    Print() << "   Writing plotfile '" << plotfilename << "'" << std::endl;

    // Write plotfile using AMReX's multi-level plotfile capability
    WriteMultiLevelPlotfile(plotfilename, finest_level+1, GetVecOfConstPtrs(mf),
                           varnames, geom, time, istep, ref_ratio);
}

void writeTecplot(const std::string& plot_file,
                  const int step,
                  const Geometry& geom,
                  const MultiFab& eta,
                  const MultiFab& u,
                  const MultiFab& v,
                  const BoxArray& ba_root,
                  const DistributionMapping& dmap_root,
                  const Real time)
{
    std::string tecfilename = Concatenate(plot_file + "_tecplot", step, 5);
    tecfilename += ".dat";

    Print() << "   Writing Tecplot file '" << tecfilename << "'" << std::endl;

    // Gather data to root processor
    const Box& domain = geom.Domain();
    const int nx = domain.length(0);
    const int ny = domain.length(1);
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    // Create MultiFabs on root processor
    MultiFab eta_root, u_root, v_root;
    eta_root.define(ba_root, dmap_root, 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    u_root.define(ba_root, dmap_root, 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    v_root.define(ba_root, dmap_root, 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

    // Gather data to root using parallel copy
    eta_root.ParallelCopy(eta, 0, 0, 1);
    u_root.ParallelCopy(u, 0, 0, 1);
    v_root.ParallelCopy(v, 0, 0, 1);

    // Only root processor writes the file
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream tecfile;
        openOutputFileOrAbort(tecfile, tecfilename, "Tecplot file");
        tecfile << std::scientific << std::setprecision(16);

        // Write Tecplot header
        tecfile << "TITLE = \"HOS Wave Surface Data\"\n";
        tecfile << "VARIABLES = \"X\", \"Y\", \"eta\", \"u\", \"v\"\n";
        tecfile << "ZONE T=\"Time=" << time << "\", I=" << nx << ", J=" << ny
                << ", F=POINT\n";

        // Write data in point format
        for (MFIter mfi(eta_root); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            Array4<Real const> eta_arr = eta_root.const_array(mfi);
            Array4<Real const> u_arr = u_root.const_array(mfi);
            Array4<Real const> v_arr = v_root.const_array(mfi);

            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    // Physical coordinates at collocation points (cell centers for spectral methods)
                    Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    Real y = prob_lo[1] + (j + 0.5) * dx[1];

                    tecfile << x << " " << y << " "
                           << eta_arr(i,j,0,0) << " "
                           << u_arr(i,j,0,0) << " "
                           << v_arr(i,j,0,0) << "\n";
                }
            }
        }

        tecfile.close();
    }

    ParallelDescriptor::Barrier();
}

void writeTecplotMultiLevelCRESCENT(const std::string& plot_file,
                                   const int step,
                                   const int finest_level,
                                   const Vector<Geometry>& geom,
                                   const Vector<MultiFab>& eta,
                                   const Vector<MultiFab>& phi,
                                   const Vector<MultiFab>& u,
                                   const Vector<MultiFab>& v,
                                   const Vector<BoxArray>& ba_root,
                                   const Vector<DistributionMapping>& dmap_root,
                                   const Real time,
                                   const NestedTecplotDiagnosticLevels& diagnostic_levels)
{
    const auto& core_surface_derivatives = diagnostic_levels.core_surface_derivatives;
    const auto& nested_output_derivatives = diagnostic_levels.nested_output_derivatives;
    const auto& chi_D_mask = diagnostic_levels.chi_D_mask;
    const auto& sigma_blend = diagnostic_levels.sigma_blend;
    const auto& eta_x_core_reference_l1 = diagnostic_levels.eta_x_core_reference_l1;
    const auto& phi_x_core_reference_l1 = diagnostic_levels.phi_x_core_reference_l1;

    std::string tecfilename = Concatenate(plot_file + "_tecplot", step, 5);
    tecfilename += ".dat";

    Print() << "   Writing multi-level Tecplot file (with CRESCENT diagnostics) '" << tecfilename << "'" << std::endl;

    // Only root processor opens the file and writes header
    std::ofstream tecfile;
    if (ParallelDescriptor::IOProcessor()) {
        openOutputFileOrAbort(tecfile, tecfilename, "CRESCENT Tecplot file");
        tecfile << std::scientific << std::setprecision(16);

        // Write Tecplot header (once for all zones)
        tecfile << "TITLE = \"HOS Wave Surface Data - Nested AMR with CRESCENT\"\n";
        tecfile << "VARIABLES = \"X\", \"Y\", \"eta\", \"phi\", \"u\", \"v\"";

        // Add CRESCENT variables for L1+
        tecfile << ", \"eta_x\", \"eta_y\", \"phi_x\", \"phi_y\", \"phi_z\", \"chi_D\", \"sigma\"";
        tecfile << ", \"eta_x_interp\", \"phi_x_interp\"";  // L0 interpolated derivatives for comparison
        tecfile << "\n";
    }

    // Write each level as a separate zone
    for (int lev = 0; lev <= finest_level; ++lev) {
        const Box& domain_box = eta[lev].boxArray().minimalBox();
        const int nx = domain_box.length(0);
        const int ny = domain_box.length(1);
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();

        // Get physical extent of this level
        const Real x_lo = prob_lo[0] + domain_box.smallEnd(0) * dx[0];
        const Real x_hi = prob_lo[0] + (domain_box.bigEnd(0) + 1) * dx[0];
        const Real y_lo = prob_lo[1] + domain_box.smallEnd(1) * dx[1];
        const Real y_hi = prob_lo[1] + (domain_box.bigEnd(1) + 1) * dx[1];

        if (ParallelDescriptor::IOProcessor()) {
            Print() << "     Level " << lev << ": " << nx << " x " << ny
                    << " cells, extent [" << x_lo << ", " << x_hi << "] x ["
                    << y_lo << ", " << y_hi << "]" << std::endl;

            // Write zone header
            tecfile << "ZONE T=\"Level " << lev << " (Time=" << time << ")\", ";
            tecfile << "I=" << nx << ", J=" << ny << ", ";
            tecfile << "F=POINT, ";
            tecfile << "SOLUTIONTIME=" << time << "\n";
        }

        // Gather data to root processor (ALL processors participate)
        MultiFab eta_root, phi_root, u_root, v_root;
        eta_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        phi_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        u_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        v_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

        eta_root.ParallelCopy(eta[lev], 0, 0, 1);
        phi_root.ParallelCopy(phi[lev], 0, 0, 1);
        u_root.ParallelCopy(u[lev], 0, 0, 1);
        v_root.ParallelCopy(v[lev], 0, 0, 1);

        // Gather derivatives (for ALL levels) and CRESCENT diagnostics (L1+ only)
        MultiFab eta_x_root, eta_y_root, phi_x_root, phi_y_root, phi_z_root;
        MultiFab chi_D_root, sigma_root;
        MultiFab eta_x_core_interp_root, phi_x_core_interp_root;

        // Choose correct derivative source:
        // Core: use standard core-path derivatives.
        // Nested L1: use the CRESCENT-blended derivatives (eta_x_RK, etc.) in Omega_M.
        const Vector<MultiFab>& eta_x_src =
            (lev == 0) ? core_surface_derivatives.eta_x : nested_output_derivatives.eta_x;
        const Vector<MultiFab>& eta_y_src =
            (lev == 0) ? core_surface_derivatives.eta_y : nested_output_derivatives.eta_y;
        const Vector<MultiFab>& phi_x_src =
            (lev == 0) ? core_surface_derivatives.phi_x : nested_output_derivatives.phi_x;
        const Vector<MultiFab>& phi_y_src =
            (lev == 0) ? core_surface_derivatives.phi_y : nested_output_derivatives.phi_y;
        const Vector<MultiFab>& phi_z_src =
            (lev == 0) ? core_surface_derivatives.phi_z : nested_output_derivatives.phi_z;

        // Check if derivatives exist for this level
        bool has_output_derivatives = (lev < (int)eta_x_src.size());
        // CRESCENT diagnostics (chi_D, sigma) only exist for L1+
        bool has_crescent_diagnostics = (lev >= 1 && lev < (int)chi_D_mask.size());
        // Interpolated derivatives only exist for L1+
        bool has_core_reference_derivatives =
            (lev >= 1 && lev < (int)eta_x_core_reference_l1.size());

        // Gather derivatives (use correct source based on level)
        if (has_output_derivatives) {
            eta_x_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
            eta_y_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
            phi_x_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
            phi_y_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
            phi_z_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

            eta_x_root.ParallelCopy(eta_x_src[lev], 0, 0, 1);
            eta_y_root.ParallelCopy(eta_y_src[lev], 0, 0, 1);
            phi_x_root.ParallelCopy(phi_x_src[lev], 0, 0, 1);
            phi_y_root.ParallelCopy(phi_y_src[lev], 0, 0, 1);
            phi_z_root.ParallelCopy(phi_z_src[lev], 0, 0, 1);
        }

        // Gather CRESCENT diagnostics (L1+ only)
        if (has_crescent_diagnostics) {
            chi_D_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
            sigma_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

            chi_D_root.ParallelCopy(chi_D_mask[lev], 0, 0, 1);
            sigma_root.ParallelCopy(sigma_blend[lev], 0, 0, 1);
        }

        // Gather interpolated derivatives (L1+ only)
        if (has_core_reference_derivatives) {
            eta_x_core_interp_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
            phi_x_core_interp_root.define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));

            eta_x_core_interp_root.ParallelCopy(eta_x_core_reference_l1[lev], 0, 0, 1);
            phi_x_core_interp_root.ParallelCopy(phi_x_core_reference_l1[lev], 0, 0, 1);
        }

        // Only root processor writes data
        if (ParallelDescriptor::IOProcessor()) {
            // Write data in point format (structured grid)
            for (MFIter mfi(eta_root); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                Array4<Real const> eta_arr = eta_root.const_array(mfi);
                Array4<Real const> phi_arr = phi_root.const_array(mfi);
                Array4<Real const> u_arr = u_root.const_array(mfi);
                Array4<Real const> v_arr = v_root.const_array(mfi);

                // Get derivative arrays (for all levels)
                Array4<Real const> eta_x_arr, eta_y_arr, phi_x_arr, phi_y_arr, phi_z_arr;
                if (has_output_derivatives) {
                    eta_x_arr = eta_x_root.const_array(mfi);
                    eta_y_arr = eta_y_root.const_array(mfi);
                    phi_x_arr = phi_x_root.const_array(mfi);
                    phi_y_arr = phi_y_root.const_array(mfi);
                    phi_z_arr = phi_z_root.const_array(mfi);
                }

                // Get CRESCENT diagnostic arrays (L1+ only)
                Array4<Real const> chi_D_arr, sigma_arr;
                if (has_crescent_diagnostics) {
                    chi_D_arr = chi_D_root.const_array(mfi);
                    sigma_arr = sigma_root.const_array(mfi);
                }

                // Get interpolated derivative arrays (L1+ only)
                Array4<Real const> eta_x_core_interp_arr, phi_x_core_interp_arr;
                if (has_core_reference_derivatives) {
                    eta_x_core_interp_arr = eta_x_core_interp_root.const_array(mfi);
                    phi_x_core_interp_arr = phi_x_core_interp_root.const_array(mfi);
                }

                // The domain box only, not the padded root box.
                Box data_box = bx & domain_box;

                for (int j = data_box.smallEnd(1); j <= data_box.bigEnd(1); ++j) {
                    for (int i = data_box.smallEnd(0); i <= data_box.bigEnd(0); ++i) {
                        // Physical coordinates at collocation points
                        Real x = prob_lo[0] + (i + 0.5) * dx[0];
                        Real y = prob_lo[1] + (j + 0.5) * dx[1];

                        // Write basic variables
                        tecfile << x << " " << y << " "
                               << eta_arr(i,j,0,0) << " "
                               << phi_arr(i,j,0,0) << " "
                               << u_arr(i,j,0,0) << " "
                               << v_arr(i,j,0,0);

                        // Write derivatives (for ALL levels)
                        if (has_output_derivatives) {
                            tecfile << " " << eta_x_arr(i,j,0,0)
                                   << " " << eta_y_arr(i,j,0,0)
                                   << " " << phi_x_arr(i,j,0,0)
                                   << " " << phi_y_arr(i,j,0,0)
                                   << " " << phi_z_arr(i,j,0,0);
                        } else {
                            // Should not happen, but write zeros as fallback
                            tecfile << " 0.0 0.0 0.0 0.0 0.0";
                        }

                        // Write CRESCENT diagnostics (L1+ only, zeros for L0)
                        if (has_crescent_diagnostics) {
                            tecfile << " " << chi_D_arr(i,j,0,0)
                                   << " " << sigma_arr(i,j,0,0);
                        } else {
                            // L0: no CRESCENT diagnostics (chi_D, sigma)
                            tecfile << " 0.0 0.0";
                        }

                        // Write interpolated derivatives (L1+ only, zeros for L0)
                        if (has_core_reference_derivatives) {
                            tecfile << " " << eta_x_core_interp_arr(i,j,0,0)
                                   << " " << phi_x_core_interp_arr(i,j,0,0);
                        } else {
                            // L0: no interpolated derivatives
                            tecfile << " 0.0 0.0";
                        }

                        tecfile << "\n";
                    }
                }
            }
        }

        ParallelDescriptor::Barrier(); // Sync after each level
    }

    if (ParallelDescriptor::IOProcessor()) {
        tecfile.close();
        Print() << "   Tecplot file written successfully with " << (finest_level+1)
                << " levels (CRESCENT diagnostics included for L1+)" << std::endl;
    }

    ParallelDescriptor::Barrier();
}

void writeTecplotMovie(const std::string& plot_file,
                       const int step,
                       const int finest_level,
                       const Vector<Geometry>& geom,
                       const Vector<MultiFab>& eta,
                       const Vector<BoxArray>& ba_root,
                       const Vector<DistributionMapping>& dmap_root,
                       const Real time,
                       const Vector<crescent::nested::domains::DomainStructure>& crescent_domains)
{
    std::string tecfilename = Concatenate(plot_file + "_movie_tecplot", step, 5);
    tecfilename += ".dat";

    Print() << "   Writing movie Tecplot file (eta only, L1=Omega_A) '" << tecfilename << "'" << std::endl;

    // Gather data to IO processor
    Vector<MultiFab> eta_root(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        eta_root[lev].define(ba_root[lev], dmap_root[lev], 1, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        eta_root[lev].ParallelCopy(eta[lev], 0, 0, 1);
    }

    // Only root processor opens the file and writes
    std::ofstream tecfile;
    if (ParallelDescriptor::IOProcessor()) {
        openOutputFileOrAbort(tecfile, tecfilename, "movie Tecplot file");
        tecfile << std::scientific << std::setprecision(16);

        // Write Tecplot header (once for all zones)
        tecfile << "TITLE = \"HOS Wave Surface Data - Movie Output (eta only, L1=Omega_A)\"\n";
        tecfile << "VARIABLES = \"X\", \"Y\", \"eta\"\n";
    }

    // Write each level as a separate zone
    for (int lev = 0; lev <= finest_level; ++lev) {
        // Determine output region
        Box output_box;
        std::string zone_label;

        if (lev == 0) {
            // L0: Output full domain
            output_box = eta[lev].boxArray().minimalBox();
            zone_label = "Level 0 (Full Domain)";
        } else {
            // L1: Output only Ω_A region
            output_box = crescent_domains[lev].box_A;
            zone_label = "Level 1 (Omega_A Region)";
        }

        const int nx = output_box.length(0);
        const int ny = output_box.length(1);
        const auto dx = geom[lev].CellSizeArray();
        const auto prob_lo = geom[lev].ProbLoArray();

        // Get physical extent of output region
        const Real x_lo = prob_lo[0] + output_box.smallEnd(0) * dx[0];
        const Real x_hi = prob_lo[0] + (output_box.bigEnd(0) + 1) * dx[0];
        const Real y_lo = prob_lo[1] + output_box.smallEnd(1) * dx[1];
        const Real y_hi = prob_lo[1] + (output_box.bigEnd(1) + 1) * dx[1];

        if (ParallelDescriptor::IOProcessor()) {
            Print() << "     Level " << lev << ": " << nx << " x " << ny
                    << " cells, extent [" << x_lo << ", " << x_hi << "] x ["
                    << y_lo << ", " << y_hi << "]" << std::endl;

            // Write ZONE header
            tecfile << "ZONE T=\"" << zone_label << " (Time="
                    << std::scientific << std::setprecision(16) << time
                    << ")\", I=" << nx << ", J=" << ny
                    << ", F=POINT, SOLUTIONTIME=" << time << "\n";
        }

        // Only root processor writes data
        if (ParallelDescriptor::IOProcessor()) {
            // Write data in point format (structured grid)
            for (MFIter mfi(eta_root[lev]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                Array4<Real const> eta_arr = eta_root[lev].const_array(mfi);

                // Get intersection with output region
                Box data_box = bx & output_box;

                if (data_box.ok()) {
                    // Loop over output region only
                    for (int j = data_box.smallEnd(1); j <= data_box.bigEnd(1); ++j) {
                        for (int i = data_box.smallEnd(0); i <= data_box.bigEnd(0); ++i) {
                            // Physical coordinates at cell centers
                            Real x = prob_lo[0] + (i + 0.5) * dx[0];
                            Real y = prob_lo[1] + (j + 0.5) * dx[1];
                            Real eta_val = eta_arr(i, j, 0);

                            // Write in scientific notation
                            tecfile << std::scientific << std::setprecision(16)
                                    << x << " " << y << " " << eta_val << "\n";
                        }
                    }
                }
            }
        }

        ParallelDescriptor::Barrier(); // Sync after each level
    }

    if (ParallelDescriptor::IOProcessor()) {
        tecfile.close();
        Print() << "   Movie Tecplot file written successfully with " << (finest_level+1)
                << " levels (eta only, L1 restricted to Omega_A)" << std::endl;
    }

    ParallelDescriptor::Barrier();
}

} // namespace crescent::io

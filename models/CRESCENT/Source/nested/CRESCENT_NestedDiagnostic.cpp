#include "CRESCENT_NestedDiagnostic.H"
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <iomanip>
#include <fstream>

namespace crescent::nested::diagnostic {

using namespace amrex;

void writeDiagnosticTecplot(
    const std::string& filename,
    const Geometry& geom,
    const MultiFab& eta_interp,
    const MultiFab& eta_periodic,
    const MultiFab& eta_x_interp,
    const MultiFab& eta_periodic_x,
    const MultiFab& eta_y_interp,
    const MultiFab& eta_periodic_y,
    const MultiFab& phi_interp,
    const MultiFab& phi_periodic,
    const MultiFab& phi_x_interp,
    const MultiFab& phi_periodic_x,
    const MultiFab& phi_x_blended,
    const MultiFab& phi_y_interp,
    const MultiFab& phi_periodic_y,
    const MultiFab& phi_y_blended,
    const MultiFab& phi_z_L1,
    const MultiFab& phi_z_L0_interp,
    const MultiFab& phi_z_blended,
    const BoxArray& ba_root,
    const DistributionMapping& dmap_root)
{
    Print() << "\n   Writing CRESCENT diagnostic Tecplot file '" << filename << "'" << std::endl;

    const Box& domain_box = phi_interp.boxArray().minimalBox();
    const int nx = domain_box.length(0);
    const int ny = domain_box.length(1);
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    // Gather data to root processor
    MultiFab eta_interp_root, eta_periodic_root;
    MultiFab eta_x_interp_root, eta_periodic_x_root, eta_y_interp_root, eta_periodic_y_root;
    MultiFab phi_interp_root, phi_periodic_root;
    MultiFab phi_x_interp_root, phi_periodic_x_root, phi_x_blended_root;
    MultiFab phi_y_interp_root, phi_periodic_y_root, phi_y_blended_root;
    MultiFab phi_z_L1_root, phi_z_L0_interp_root, phi_z_blended_root;

    eta_interp_root.define(ba_root, dmap_root, 1, 0);
    eta_periodic_root.define(ba_root, dmap_root, 1, 0);
    eta_x_interp_root.define(ba_root, dmap_root, 1, 0);
    eta_periodic_x_root.define(ba_root, dmap_root, 1, 0);
    eta_y_interp_root.define(ba_root, dmap_root, 1, 0);
    eta_periodic_y_root.define(ba_root, dmap_root, 1, 0);
    phi_interp_root.define(ba_root, dmap_root, 1, 0);
    phi_periodic_root.define(ba_root, dmap_root, 1, 0);
    phi_x_interp_root.define(ba_root, dmap_root, 1, 0);
    phi_periodic_x_root.define(ba_root, dmap_root, 1, 0);
    phi_x_blended_root.define(ba_root, dmap_root, 1, 0);
    phi_y_interp_root.define(ba_root, dmap_root, 1, 0);
    phi_periodic_y_root.define(ba_root, dmap_root, 1, 0);
    phi_y_blended_root.define(ba_root, dmap_root, 1, 0);
    phi_z_L1_root.define(ba_root, dmap_root, 1, 0);
    phi_z_L0_interp_root.define(ba_root, dmap_root, 1, 0);
    phi_z_blended_root.define(ba_root, dmap_root, 1, 0);

    eta_interp_root.ParallelCopy(eta_interp, 0, 0, 1);
    eta_periodic_root.ParallelCopy(eta_periodic, 0, 0, 1);
    eta_x_interp_root.ParallelCopy(eta_x_interp, 0, 0, 1);
    eta_periodic_x_root.ParallelCopy(eta_periodic_x, 0, 0, 1);
    eta_y_interp_root.ParallelCopy(eta_y_interp, 0, 0, 1);
    eta_periodic_y_root.ParallelCopy(eta_periodic_y, 0, 0, 1);
    phi_interp_root.ParallelCopy(phi_interp, 0, 0, 1);
    phi_periodic_root.ParallelCopy(phi_periodic, 0, 0, 1);
    phi_x_interp_root.ParallelCopy(phi_x_interp, 0, 0, 1);
    phi_periodic_x_root.ParallelCopy(phi_periodic_x, 0, 0, 1);
    phi_x_blended_root.ParallelCopy(phi_x_blended, 0, 0, 1);
    phi_y_interp_root.ParallelCopy(phi_y_interp, 0, 0, 1);
    phi_periodic_y_root.ParallelCopy(phi_periodic_y, 0, 0, 1);
    phi_y_blended_root.ParallelCopy(phi_y_blended, 0, 0, 1);
    phi_z_L1_root.ParallelCopy(phi_z_L1, 0, 0, 1);
    phi_z_L0_interp_root.ParallelCopy(phi_z_L0_interp, 0, 0, 1);
    phi_z_blended_root.ParallelCopy(phi_z_blended, 0, 0, 1);

    // Only root processor writes
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream tecfile(filename);
        if (!tecfile.is_open()) {
            amrex::Abort("Failed to open nested initialization derivative diagnostic for writing");
        }
        tecfile << std::scientific << std::setprecision(16);

        tecfile << "TITLE = \"CRESCENT Derivative Pipeline Diagnostic (C³-Repair)\"\n";

        tecfile << "VARIABLES = \"X\", \"Y\", ";
        tecfile << "\"eta_interp\", ";
        tecfile << "\"eta_periodic\", ";
        tecfile << "\"eta_x_interp\", ";
        tecfile << "\"eta_periodic_x\", ";
        tecfile << "\"eta_y_interp\", ";
        tecfile << "\"eta_periodic_y\", ";
        tecfile << "\"phi_interp\", ";
        tecfile << "\"phi_periodic\", ";
        tecfile << "\"phi_x_interp\", ";
        tecfile << "\"phi_periodic_x\", ";
        tecfile << "\"phi_x_blended\", ";
        tecfile << "\"phi_y_interp\", ";
        tecfile << "\"phi_periodic_y\", ";
        tecfile << "\"phi_y_blended\", ";
        tecfile << "\"phi_z_L1\", \"phi_z_L0_interp\", \"phi_z_blended\"\n";

        // Write zone header
        tecfile << "ZONE T=\"L1 Initialization\", ";
        tecfile << "I=" << nx << ", J=" << ny << ", ";
        tecfile << "F=POINT\n";

        // Write data
        for (MFIter mfi(phi_interp_root); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            Box data_box = bx & domain_box;

            auto eta_interp_arr = eta_interp_root.const_array(mfi);
            auto eta_periodic_arr = eta_periodic_root.const_array(mfi);
            auto eta_x_interp_arr = eta_x_interp_root.const_array(mfi);
            auto eta_periodic_x_arr = eta_periodic_x_root.const_array(mfi);
            auto eta_y_interp_arr = eta_y_interp_root.const_array(mfi);
            auto eta_periodic_y_arr = eta_periodic_y_root.const_array(mfi);
            auto phi_interp_arr = phi_interp_root.const_array(mfi);
            auto phi_periodic_arr = phi_periodic_root.const_array(mfi);
            auto phi_x_interp_arr = phi_x_interp_root.const_array(mfi);
            auto phi_periodic_x_arr = phi_periodic_x_root.const_array(mfi);
            auto phi_x_blended_arr = phi_x_blended_root.const_array(mfi);
            auto phi_y_interp_arr = phi_y_interp_root.const_array(mfi);
            auto phi_periodic_y_arr = phi_periodic_y_root.const_array(mfi);
            auto phi_y_blended_arr = phi_y_blended_root.const_array(mfi);
            auto phi_z_L1_arr = phi_z_L1_root.const_array(mfi);
            auto phi_z_L0_interp_arr = phi_z_L0_interp_root.const_array(mfi);
            auto phi_z_blended_arr = phi_z_blended_root.const_array(mfi);

            for (int j = data_box.smallEnd(1); j <= data_box.bigEnd(1); ++j) {
                for (int i = data_box.smallEnd(0); i <= data_box.bigEnd(0); ++i) {
                    Real x = prob_lo[0] + (i + 0.5) * dx[0];
                    Real y = prob_lo[1] + (j + 0.5) * dx[1];

                    tecfile << x << " " << y << " "
                           << eta_interp_arr(i,j,0) << " "
                           << eta_periodic_arr(i,j,0) << " "
                           << eta_x_interp_arr(i,j,0) << " "
                           << eta_periodic_x_arr(i,j,0) << " "
                           << eta_y_interp_arr(i,j,0) << " "
                           << eta_periodic_y_arr(i,j,0) << " "
                           << phi_interp_arr(i,j,0) << " "
                           << phi_periodic_arr(i,j,0) << " "
                           << phi_x_interp_arr(i,j,0) << " "
                           << phi_periodic_x_arr(i,j,0) << " "
                           << phi_x_blended_arr(i,j,0) << " "
                           << phi_y_interp_arr(i,j,0) << " "
                           << phi_periodic_y_arr(i,j,0) << " "
                           << phi_y_blended_arr(i,j,0) << " "
                           << phi_z_L1_arr(i,j,0) << " "
                           << phi_z_L0_interp_arr(i,j,0) << " "
                           << phi_z_blended_arr(i,j,0) << "\n";
                }
            }
        }

        tecfile.close();
        Print() << "   CRESCENT diagnostic file written successfully" << std::endl;
    }

    ParallelDescriptor::Barrier();
}

} // namespace crescent::nested::diagnostic

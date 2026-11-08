#include "CRESCENT_NestedEnergy.H"
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <fstream>
#include <iomanip>

namespace crescent::nested::energy {

using namespace amrex;

Real computeTotalEnergy(
    const MultiFab& eta,
    const MultiFab& psi,
    const MultiFab& eta_x,
    const MultiFab& eta_y,
    const MultiFab& psi_x,
    const MultiFab& psi_y,
    const MultiFab& w,
    const Geometry& geom,
    Real gravity)
{
    const Real* dx = geom.CellSize();
    Real dA = dx[0] * dx[1];  // Cell area

    // Zakharov Hamiltonian: E = ∫∫ (g·η²/2 + ψ·G(η)ψ/2) dA.
    Real energy_local = 0.0;

    for (MFIter mfi(eta); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto eta_arr = eta.const_array(mfi);
        auto psi_arr = psi.const_array(mfi);
        auto eta_x_arr = eta_x.const_array(mfi);
        auto eta_y_arr = eta_y.const_array(mfi);
        auto psi_x_arr = psi_x.const_array(mfi);
        auto psi_y_arr = psi_y.const_array(mfi);
        auto w_arr = w.const_array(mfi);

        // Use ReduceOps for proper reduction
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename ReduceData<Real>::Type;

        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            Real potential = 0.5 * gravity * eta_arr(i,j,k) * eta_arr(i,j,k);
            const Real metric = 1.0
                + eta_x_arr(i,j,k) * eta_x_arr(i,j,k)
                + eta_y_arr(i,j,k) * eta_y_arr(i,j,k);
            const Real dno = metric * w_arr(i,j,k)
                - eta_x_arr(i,j,k) * psi_x_arr(i,j,k)
                - eta_y_arr(i,j,k) * psi_y_arr(i,j,k);
            Real kinetic = 0.5 * psi_arr(i,j,k) * dno;
            return (potential + kinetic) * dA;
        });

        Real box_energy = amrex::get<0>(reduce_data.value());
        energy_local += box_energy;
    }

    // Sum across MPI ranks
    Real energy_total = energy_local;
    ParallelDescriptor::ReduceRealSum(energy_total);

    return energy_total;
}

Real computeTotalEnergyInBox(
    const MultiFab& eta,
    const MultiFab& psi,
    const MultiFab& eta_x,
    const MultiFab& eta_y,
    const MultiFab& psi_x,
    const MultiFab& psi_y,
    const MultiFab& w,
    const Geometry& geom,
    const Box& region,
    Real gravity)
{
    const Real* dx = geom.CellSize();
    Real dA = dx[0] * dx[1];

    Real energy_local = 0.0;

    for (MFIter mfi(eta); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox() & region;
        if (bx.isEmpty()) {
            continue;
        }
        auto eta_arr = eta.const_array(mfi);
        auto psi_arr = psi.const_array(mfi);
        auto eta_x_arr = eta_x.const_array(mfi);
        auto eta_y_arr = eta_y.const_array(mfi);
        auto psi_x_arr = psi_x.const_array(mfi);
        auto psi_y_arr = psi_y.const_array(mfi);
        auto w_arr = w.const_array(mfi);

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename ReduceData<Real>::Type;

        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
            Real potential = 0.5 * gravity * eta_arr(i,j,k) * eta_arr(i,j,k);
            const Real metric = 1.0
                + eta_x_arr(i,j,k) * eta_x_arr(i,j,k)
                + eta_y_arr(i,j,k) * eta_y_arr(i,j,k);
            const Real dno = metric * w_arr(i,j,k)
                - eta_x_arr(i,j,k) * psi_x_arr(i,j,k)
                - eta_y_arr(i,j,k) * psi_y_arr(i,j,k);
            Real kinetic = 0.5 * psi_arr(i,j,k) * dno;
            return (potential + kinetic) * dA;
        });

        energy_local += amrex::get<0>(reduce_data.value());
    }

    Real energy_total = energy_local;
    ParallelDescriptor::ReduceRealSum(energy_total);

    return energy_total;
}

void writeEnergyToCSV(
    const std::string& filename,
    int istep,
    Real time,
    Real energy_L0,
    Real energy_L1,
    Real energy_initial_total,
    Real energy_L1_omega_a,
    Real energy_L0_under_omega_a,
    Real energy_composite)
{
    // Only write from rank 0
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream ofs;

        // Check if file exists
        std::ifstream test(filename);
        bool file_exists = test.good();
        test.close();

        // Open in append mode
        ofs.open(filename, std::ios::app);

        // Write header if new file
        if (!file_exists) {
            ofs << "timestep,time,energy_L0,energy_L1,energy_total,energy_change_percent,"
                   "energy_L1_omegaA,energy_L0_under_omegaA,energy_composite\n";
        }

        Real energy_total = energy_L0 + energy_L1;

        const Real energy_change_pct = (energy_initial_total > 0.0)
            ? 100.0 * (energy_total - energy_initial_total) / energy_initial_total
            : 0.0;

        ofs << std::setprecision(6) << std::scientific
            << istep << ","
            << time << ","
            << energy_L0 << ","
            << energy_L1 << ","
            << energy_total << ","
            << std::fixed << std::setprecision(4) << energy_change_pct << ","
            << std::setprecision(6) << std::scientific
            << energy_L1_omega_a << ","
            << energy_L0_under_omega_a << ","
            << energy_composite << "\n";

        ofs.close();
    }
}

} // namespace crescent::nested::energy

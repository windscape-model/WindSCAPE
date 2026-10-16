#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <iostream>

#include "CRESCENT_NestedEnergy.H"

namespace {

bool nearlyEqual(amrex::Real lhs, amrex::Real rhs)
{
    constexpr amrex::Real tolerance = 1.0e-13;
    return std::abs(lhs - rhs) <= tolerance * (1.0 + std::abs(rhs));
}

int runTest()
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(1, 1, 0));
    const RealBox real_box({0.0, 0.0, 0.0}, {2.0, 2.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 0};
    const Geometry geom(domain, &real_box, 0, periodic);
    const BoxArray grids(domain);
    const DistributionMapping mapping(grids);

    MultiFab eta(grids, mapping, 1, 0);
    MultiFab psi(grids, mapping, 1, 0);
    MultiFab eta_x(grids, mapping, 1, 0);
    MultiFab eta_y(grids, mapping, 1, 0);
    MultiFab psi_x(grids, mapping, 1, 0);
    MultiFab psi_y(grids, mapping, 1, 0);
    MultiFab w(grids, mapping, 1, 0);

    eta.setVal(2.0);
    psi.setVal(3.0);
    eta_x.setVal(0.2);
    eta_y.setVal(0.1);
    psi_x.setVal(0.4);
    psi_y.setVal(0.5);
    w.setVal(0.6);

    constexpr Real gravity = 9.81;
    // G(eta)psi = (1 + 0.2^2 + 0.1^2) 0.6 - 0.2*0.4 - 0.1*0.5 = 0.5.
    constexpr Real energy_per_cell = 0.5 * gravity * 4.0 + 0.5 * 3.0 * 0.5;
    constexpr Real expected_total = 4.0 * energy_per_cell;

    const Real total = crescent::nested::energy::computeTotalEnergy(
        eta, psi, eta_x, eta_y, psi_x, psi_y, w, geom, gravity);
    const Box one_cell(IntVect(0, 0, 0), IntVect(0, 0, 0));
    const Real restricted = crescent::nested::energy::computeTotalEnergyInBox(
        eta, psi, eta_x, eta_y, psi_x, psi_y, w, geom, one_cell, gravity);

    if (!nearlyEqual(total, expected_total) ||
        !nearlyEqual(restricted, energy_per_cell)) {
        if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "Zakharov energy mismatch: total=" << total
                      << " expected=" << expected_total
                      << ", restricted=" << restricted
                      << " expected=" << energy_per_cell << '\n';
        }
        return 1;
    }
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = runTest();
    amrex::Finalize();
    return result;
}

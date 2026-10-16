#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <iostream>

#include "CRESCENT_Interpolation.H"

namespace {

bool nearlyEqual(amrex::Real lhs, amrex::Real rhs)
{
    constexpr amrex::Real tolerance = 1.0e-13;
    return std::abs(lhs - rhs) <= tolerance * (1.0 + std::abs(rhs));
}

bool testRefinementRatio(int ratio)
{
    using namespace amrex;

    const Box coarse_domain(IntVect(0, 0, 0), IntVect(7, 7, 0));
    const Box fine_domain(
        IntVect(0, 0, 0),
        IntVect(8 * ratio - 1, 8 * ratio - 1, 0));
    const Box fine_test_box(
        IntVect(ratio, ratio, 0),
        IntVect(7 * ratio - 1, 7 * ratio - 1, 0));
    const RealBox real_box({0.0, 0.0, 0.0}, {8.0, 8.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 0};
    const Geometry coarse_geom(coarse_domain, &real_box, 0, periodic);
    const Geometry fine_geom(fine_domain, &real_box, 0, periodic);
    const BoxArray coarse_grids(coarse_domain);
    const BoxArray fine_grids(fine_test_box);
    const DistributionMapping coarse_mapping(coarse_grids);
    const DistributionMapping fine_mapping(fine_grids);
    const IntVect ref_ratio(ratio, ratio, 1);

    MultiFab coarse(coarse_grids, coarse_mapping, 1, 0);
    MultiFab corrected(fine_grids, fine_mapping, 1, 0);
    MultiFab exact(fine_grids, fine_mapping, 1, 0);
    MultiFab error(fine_grids, fine_mapping, 1, 0);
    MultiFab coarse_patch(
        amrex::coarsen(fine_grids, ref_ratio), fine_mapping, 1, 1);

    for (MFIter mfi(coarse); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto field = coarse.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = Real(i) + Real(0.5);
            const Real y = Real(j) + Real(0.5);
            field(i, j, k) = Real(2.0) * x - Real(3.0) * y + Real(5.0);
        });
    }
    for (MFIter mfi(exact); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto field = exact.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = (Real(i) + Real(0.5)) / Real(ratio);
            const Real y = (Real(j) + Real(0.5)) / Real(ratio);
            field(i, j, k) = Real(2.0) * x - Real(3.0) * y + Real(5.0);
        });
    }

    crescent::interpolation::interpolateScalarParentToChild(
        coarse,
        corrected,
        coarse_geom,
        fine_geom,
        ref_ratio,
        coarse_patch,
        1);

    // The cell-centred stencil reproduces an affine field exactly.
    MultiFab::Copy(error, corrected, 0, 0, 1, 0);
    MultiFab::Subtract(error, exact, 0, 0, 1, 0);
    const Real corrected_error = error.norm0();
    if (!nearlyEqual(corrected_error, 0.0)) {
        if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "bilinear interpolation error at ratio " << ratio
                      << ": " << corrected_error << '\n';
        }
        return false;
    }

    return true;
}

int runTest()
{
    return testRefinementRatio(2) && testRefinementRatio(4) ? 0 : 1;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = runTest();
    amrex::Finalize();
    return result;
}

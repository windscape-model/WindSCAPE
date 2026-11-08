#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <iostream>
#include <memory>

#include "CRESCENT_Dynamics.H"

namespace {

bool nearlyEqual(amrex::Real lhs, amrex::Real rhs)
{
    constexpr amrex::Real tolerance = 1.0e-13;
    return std::abs(lhs - rhs) <= tolerance * (1.0 + std::abs(rhs));
}

int runTest()
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(7, 7, 0));
    const RealBox real_box({0.0, 0.0, 0.0}, {8.0, 8.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 0};
    const Geometry geom(domain, &real_box, 0, periodic);
    const BoxArray grids(domain);
    const DistributionMapping mapping(grids);

    MultiFab eta(grids, mapping, 1, 0);
    MultiFab phi(grids, mapping, 1, 0);
    MultiFab eta_x(grids, mapping, 1, 0);
    MultiFab eta_y(grids, mapping, 1, 0);
    MultiFab phi_x(grids, mapping, 1, 0);
    MultiFab phi_y(grids, mapping, 1, 0);
    MultiFab phi_z(grids, mapping, 1, 0);
    MultiFab pressure(grids, mapping, 1, 0);
    MultiFab rhs_eta(grids, mapping, 4, 0);
    MultiFab rhs_phi(grids, mapping, 4, 0);
    MultiFab metric_scratch(grids, mapping, 1, 0);
    MultiFab phi_z_sq_scratch(grids, mapping, 1, 0);
    MultiFab dealiasing_scratch(grids, mapping, 1, 0);

    constexpr Real eta_value = 2.0;
    constexpr Real phi_z_value = 3.0;
    constexpr Real pressure_value = 0.5;
    constexpr Real gravity = 9.81;
    constexpr Real sentinel = -777.0;

    eta.setVal(eta_value);
    phi.setVal(0.0);
    eta_x.setVal(0.0);
    eta_y.setVal(0.0);
    phi_x.setVal(0.0);
    phi_y.setVal(0.0);
    phi_z.setVal(phi_z_value);
    pressure.setVal(pressure_value);
    rhs_eta.setVal(sentinel);
    rhs_phi.setVal(sentinel);

    std::unique_ptr<FFT::R2C<Real>> unused_fft;
    for (int stage = 0; stage < 4; ++stage) {
        crescent::dynamics::computeSurfaceRHS(
            rhs_eta,
            rhs_phi,
            eta,
            eta_x,
            eta_y,
            phi_x,
            phi_y,
            phi_z,
            pressure,
            metric_scratch,
            phi_z_sq_scratch,
            unused_fft,
            geom.data(),
            1.0 / gravity,
            1,
            stage,
            dealiasing_scratch);

        for (int component = 0; component < 4; ++component) {
            const Real expected_eta = component <= stage ? phi_z_value : sentinel;
            const Real expected_phi = component <= stage
                ? -gravity * eta_value - pressure_value
                : sentinel;
            if (!nearlyEqual(rhs_eta.min(component), expected_eta) ||
                !nearlyEqual(rhs_eta.max(component), expected_eta) ||
                !nearlyEqual(rhs_phi.min(component), expected_phi) ||
                !nearlyEqual(rhs_phi.max(component), expected_phi)) {
                if (ParallelDescriptor::IOProcessor()) {
                    std::cerr << "linear RHS mismatch after stage " << stage
                              << " in component " << component << '\n';
                }
                return 1;
            }
        }
    }

    // Nonlinear RHS pressure must be applied after the spectral projection on
    // the corrected path. A checkerboard pressure is deliberately outside the
    // retained two-thirds band, so projecting it would erase this signal.
    eta.setVal(0.0);
    phi_z.setVal(0.0);
    rhs_eta.setVal(sentinel);
    rhs_phi.setVal(sentinel);
    for (MFIter mfi(pressure); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto pa = pressure.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pa(i, j, k) = ((i + j) & 1) == 0 ? Real(1.0) : Real(-1.0);
        });
    }
    auto fft = std::make_unique<FFT::R2C<Real>>(domain);
    crescent::dynamics::computeSurfaceRHS(
        rhs_eta,
        rhs_phi,
        eta,
        eta_x,
        eta_y,
        phi_x,
        phi_y,
        phi_z,
        pressure,
        metric_scratch,
        phi_z_sq_scratch,
        fft,
        geom.data(),
        1.0 / gravity,
        2,
        0,
        dealiasing_scratch);

    MultiFab pressure_error(grids, mapping, 1, 0);
    MultiFab::Copy(pressure_error, rhs_phi, 0, 0, 1, 0);
    MultiFab::Add(pressure_error, pressure, 0, 0, 1, 0);
    if (!nearlyEqual(rhs_eta.norm0(0), 0.0) ||
        !nearlyEqual(pressure_error.norm0(0), 0.0)) {
        if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "pressure-outside-projection RHS mismatch\n";
        }
        return 1;
    }

    // The adjustment must ramp only the nonlinear remainder. With a constant
    // horizontal velocity, F_psi = -2 - p and F_psi,lin = -p, so the exact
    // adjusted result is -2 W(t) - p even for an unresolved pressure field.
    constexpr Real stage_time = 0.5;
    constexpr Real adjustment_time = 1.0;
    const Real adjustment_weight =
        crescent::dynamics::computeAdjustmentRampWeight(
            stage_time, adjustment_time);
    if (!nearlyEqual(
            adjustment_weight,
            Real(1.0) - std::exp(-Real(0.0625))) ||
        !nearlyEqual(
            crescent::dynamics::computeAdjustmentRampWeight(0.0, 1.0),
            0.0) ||
        !nearlyEqual(
            crescent::dynamics::computeAdjustmentRampWeight(2.0, 1.0),
            1.0)) {
        if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "adjustment ramp weight mismatch\n";
        }
        return 1;
    }

    phi_x.setVal(2.0);
    rhs_eta.setVal(sentinel);
    rhs_phi.setVal(sentinel);
    crescent::dynamics::computeSurfaceRHS(
        rhs_eta,
        rhs_phi,
        eta,
        eta_x,
        eta_y,
        phi_x,
        phi_y,
        phi_z,
        pressure,
        metric_scratch,
        phi_z_sq_scratch,
        fft,
        geom.data(),
        1.0 / gravity,
        2,
        0,
        dealiasing_scratch);
    crescent::dynamics::applyAdjustmentRampToStageRHS(
        rhs_eta,
        rhs_phi,
        eta,
        phi,
        pressure,
        fft,
        geom.data(),
        stage_time,
        adjustment_time,
        0,
        dealiasing_scratch);

    MultiFab adjusted_error_field(grids, mapping, 1, 0);
    for (MFIter mfi(rhs_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto values = rhs_phi.const_array(mfi);
        auto pa = pressure.const_array(mfi);
        auto error = adjusted_error_field.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real expected = -Real(2.0) * adjustment_weight
                                  - pa(i, j, k);
            error(i, j, k) = values(i, j, k, 0) - expected;
        });
    }
    const Real adjusted_error = adjusted_error_field.norm0(0);
    if (!nearlyEqual(adjusted_error, 0.0)) {
        if (ParallelDescriptor::IOProcessor()) {
            std::cerr << "adjusted nonlinear RHS mismatch: "
                      << adjusted_error << '\n';
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

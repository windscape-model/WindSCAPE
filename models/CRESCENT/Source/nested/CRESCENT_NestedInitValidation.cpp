#include "CRESCENT_NestedInitValidation.H"

// Optional nested-L1 initialization diagnostics and validation.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>

#include "CRESCENT_NestedC3Dispatch.H"
#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_NestedDerivativeExecution.H"
#include "CRESCENT_NestedDiagnostic.H"
#include "CRESCENT_NestedTypes.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_Interpolation.H"

using namespace amrex;

namespace {

using NestedL1BoundaryDiagnostic = crescent::nested::BoundaryDiagnostic;
using NestedL1InitializationValidationScratch =
    crescent::nested::InitializationValidationScratch;
using InitializationValidationContext =
    crescent::nested::validation::InitializationValidationContext;
using LevelWaveSpectralSettings = crescent::shared::LevelWaveSpectralSettings;
using ApplyNestedValidationRepairFn =
    std::function<void(amrex::MultiFab&, amrex::MultiFab&, int)>;

struct NestedInitValidationContext
{
    int perturbation_order;
    int verbosity;
    int c3_use_interface;
    int linear_derivative_spectral_filter_mode;
    int derivative_interpolation_order;
    LevelWaveSpectralSettings spectral_settings;
    amrex::IntVect refine_ratio;

    const amrex::BoxArray& grids;
    const amrex::DistributionMapping& dmap;
    const amrex::Geometry& geom;
    const amrex::Geometry& core_geom;
    const amrex::Box& box_S;
    const amrex::BoxArray& root_ba;
    const amrex::DistributionMapping& root_dmap;
    const crescent::nested::domains::DomainStructure& domain;
    int c3_buffer_width;
    crescent::shared::LevelWaveOperatorView level_wave_operator;

    amrex::MultiFab& eta;
    amrex::MultiFab& phi;
    const amrex::MultiFab& phi_x;
    const amrex::MultiFab& phi_y;
    amrex::MultiFab& wavenumber_pow;
    amrex::MultiFab& eta_pow;
    amrex::MultiFab& phi_expansion;

    const amrex::MultiFab& phi_z_core;
    const amrex::MultiFab& sigma_blend;
    std::string diagnostic_output_file;

    ApplyNestedValidationRepairFn apply_nested_c3_repair_pair;
};

NestedL1BoundaryDiagnostic collectNestedL1BoundaryDiagnostic(
    const amrex::BoxArray& level_grids,
    const amrex::Geometry& level_geom,
    const amrex::MultiFab& phi_field)
{
    const Box& global_box = level_grids.minimalBox();

    NestedL1BoundaryDiagnostic diagnostic;
    diagnostic.i_lo = global_box.smallEnd(0);
    diagnostic.i_hi = global_box.bigEnd(0);
    diagnostic.j_mid = (global_box.smallEnd(1) + global_box.bigEnd(1)) / 2;

    const Real* prob_lo = level_geom.ProbLo();
    const Real* dx_arr = level_geom.CellSize();
    diagnostic.x_left = prob_lo[0] + (diagnostic.i_lo + 0.5) * dx_arr[0];
    diagnostic.x_right = prob_lo[0] + (diagnostic.i_hi + 0.5) * dx_arr[0];
    diagnostic.y_mid = prob_lo[1] + (diagnostic.j_mid + 0.5) * dx_arr[1];

    for (MFIter mfi(phi_field); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& phi_arr = phi_field.const_array(mfi);
        if (bx.contains(IntVect(diagnostic.i_lo, diagnostic.j_mid, 0))) {
            diagnostic.phi_left = phi_arr(diagnostic.i_lo, diagnostic.j_mid, 0);
        }
        if (bx.contains(IntVect(diagnostic.i_hi, diagnostic.j_mid, 0))) {
            diagnostic.phi_right = phi_arr(diagnostic.i_hi, diagnostic.j_mid, 0);
        }
    }

    constexpr Real k = 0.157;
    constexpr Real omega = 1.24;
    constexpr Real A = 0.636;
    const Real phi_amp = A * omega / k;
    diagnostic.phi_left_theory = phi_amp * std::sin(k * diagnostic.x_left);
    diagnostic.phi_right_theory = phi_amp * std::sin(k * diagnostic.x_right);

    return diagnostic;
}

void printNestedL1BoundaryDiagnosticBefore(
    const NestedL1BoundaryDiagnostic& diagnostic)
{
    amrex::Print() << "\n=== BOUNDARY DIAGNOSTIC (Before C3-repair) ===" << std::endl;
    amrex::Print() << "  L1 box: i ∈ [" << diagnostic.i_lo << ", " << diagnostic.i_hi << "]" << std::endl;
    amrex::Print() << "  Physical: x ∈ [" << diagnostic.x_left << ", " << diagnostic.x_right
                   << "] m" << std::endl;
    amrex::Print() << "  Checking j=" << diagnostic.j_mid << " (y=" << diagnostic.y_mid
                   << " m)" << std::endl;
    amrex::Print() << "  φ(i=" << diagnostic.i_lo << ", x=" << diagnostic.x_left << "m) = "
                   << diagnostic.phi_left << " m²/s (theory: " << diagnostic.phi_left_theory
                   << ")" << std::endl;
    amrex::Print() << "  φ(i=" << diagnostic.i_hi << ", x=" << diagnostic.x_right << "m) = "
                   << diagnostic.phi_right << " m²/s (theory: " << diagnostic.phi_right_theory
                   << ")" << std::endl;
    amrex::Print() << "  Δφ_observed = " << (diagnostic.phi_right - diagnostic.phi_left)
                   << " m²/s" << std::endl;
    amrex::Print() << "  Δφ_theory   = "
                   << (diagnostic.phi_right_theory - diagnostic.phi_left_theory)
                   << " m²/s" << std::endl;
    amrex::Print() << "  Discrepancy: "
                   << std::abs(diagnostic.phi_right - diagnostic.phi_left)
                          / std::abs(diagnostic.phi_right_theory - diagnostic.phi_left_theory)
                   << "×" << std::endl;
}

void printNestedL1BoundaryDiagnosticAfter(
    const NestedL1BoundaryDiagnostic& before,
    const NestedL1BoundaryDiagnostic& after)
{
    amrex::Print() << "\n=== BOUNDARY DIAGNOSTIC (After C3-repair) ===" << std::endl;
    amrex::Print() << "  φ(i=" << after.i_lo << ", x=" << after.x_left << "m) = "
                   << after.phi_left << " m²/s" << std::endl;
    amrex::Print() << "  φ(i=" << after.i_hi << ", x=" << after.x_right << "m) = "
                   << after.phi_right << " m²/s" << std::endl;
    amrex::Print() << "  Δφ_after = " << (after.phi_right - after.phi_left) << " m²/s" << std::endl;
    amrex::Print() << "  Reduction: " << std::abs(before.phi_right - before.phi_left)
                   << " → " << std::abs(after.phi_right - after.phi_left) << " m²/s" << std::endl;

    const Real reduction_factor = std::abs(before.phi_right - before.phi_left)
                                  / std::max(std::abs(after.phi_right - after.phi_left), 1.0e-10);
    amrex::Print() << "  Jump reduced by: " << reduction_factor << "×" << std::endl;

    if (std::abs(after.phi_right - after.phi_left) < 0.1) {
        amrex::Print() << "  ✓ C3-repair successful (jump < 0.1 m²/s)" << std::endl;
    } else {
        amrex::Print() << "  ✗ C3-repair failed (jump still large)" << std::endl;
    }
    amrex::Print() << "==========================================\n" << std::endl;
}

void writeNestedL1PeriodicitySlice(
    const amrex::BoxArray& level_grids,
    const amrex::Geometry& level_geom,
    const amrex::MultiFab& eta_after,
    const amrex::MultiFab& phi_after,
    const amrex::MultiFab& eta_before,
    const amrex::MultiFab& phi_before)
{
    amrex::Print() << "\n=== Extracting 1D slice for periodicity check ===" << std::endl;
    const int nx_nested = level_grids[0].length(0);
    const int ny_nested = level_grids[0].length(1);

    std::ofstream slice_file("c3_slice_validation.dat");
    slice_file << "# 1D slice through middle of L1 domain" << std::endl;
    slice_file << "# i_global  x_coord  eta_before  eta_after  phi_before  phi_after" << std::endl;

    const auto dx = level_geom.CellSizeArray();
    const auto prob_lo = level_geom.ProbLoArray();
    const int j_mid = level_grids[0].smallEnd(1) + ny_nested / 2;

    for (MFIter mfi(eta_after); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (bx.smallEnd(1) > j_mid || bx.bigEnd(1) < j_mid) {
            continue;
        }

        auto eta_bef = eta_before.const_array(mfi);
        auto eta_aft = eta_after.const_array(mfi);
        auto phi_bef = phi_before.const_array(mfi);
        auto phi_aft = phi_after.const_array(mfi);

        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            const Real x = prob_lo[0] + (i + 0.5) * dx[0];

            slice_file << i << " " << x << " "
                       << eta_bef(i, j_mid, 0) << " "
                       << eta_aft(i, j_mid, 0) << " "
                       << phi_bef(i, j_mid, 0) << " "
                       << phi_aft(i, j_mid, 0) << std::endl;
        }
    }
    slice_file.close();
    amrex::Print() << "  Wrote 1D slice to 'c3_slice_validation.dat'" << std::endl;
    amrex::Print() << "  L1 extent: nx=" << nx_nested << " x ny=" << ny_nested
                   << ", j_mid=" << j_mid << std::endl;
}

void logNestedL1InterpolatedDerivativeErrors(
    const amrex::BoxArray& level_grids,
    const amrex::DistributionMapping& level_dmap,
    const NestedL1InitializationValidationScratch& scratch,
    const amrex::MultiFab& eta_x_CRESCENT,
    const amrex::MultiFab& eta_y_CRESCENT,
    const amrex::MultiFab& phi_x_CRESCENT,
    const amrex::MultiFab& phi_y_CRESCENT)
{
    amrex::Print() << "\nL0 interpolated derivatives:" << std::endl;
    amrex::Print() << "  eta_x range: [" << scratch.eta_x_core_interp.min(0) << ", "
                   << scratch.eta_x_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "  eta_y range: [" << scratch.eta_y_core_interp.min(0) << ", "
                   << scratch.eta_y_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_x range: [" << scratch.phi_x_core_interp.min(0) << ", "
                   << scratch.phi_x_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_y range: [" << scratch.phi_y_core_interp.min(0) << ", "
                   << scratch.phi_y_core_interp.max(0) << "]" << std::endl;

    MultiFab eta_x_error(level_grids, level_dmap, 1, 0);
    MultiFab phi_x_error(level_grids, level_dmap, 1, 0);
    MultiFab eta_y_error(level_grids, level_dmap, 1, 0);
    MultiFab phi_y_error(level_grids, level_dmap, 1, 0);

    MultiFab::Copy(eta_x_error, eta_x_CRESCENT, 0, 0, 1, 0);
    MultiFab::Copy(phi_x_error, phi_x_CRESCENT, 0, 0, 1, 0);
    MultiFab::Copy(eta_y_error, eta_y_CRESCENT, 0, 0, 1, 0);
    MultiFab::Copy(phi_y_error, phi_y_CRESCENT, 0, 0, 1, 0);

    MultiFab::Subtract(eta_x_error, scratch.eta_x_core_interp, 0, 0, 1, 0);
    MultiFab::Subtract(phi_x_error, scratch.phi_x_core_interp, 0, 0, 1, 0);
    MultiFab::Subtract(eta_y_error, scratch.eta_y_core_interp, 0, 0, 1, 0);
    MultiFab::Subtract(phi_y_error, scratch.phi_y_core_interp, 0, 0, 1, 0);

    amrex::Print() << "\nErrors (CRESCENT - L0_interp):" << std::endl;
    amrex::Print() << "  eta_x error range: [" << eta_x_error.min(0) << ", "
                   << eta_x_error.max(0) << "]" << std::endl;
    amrex::Print() << "  eta_y error range: [" << eta_y_error.min(0) << ", "
                   << eta_y_error.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_x error range: [" << phi_x_error.min(0) << ", "
                   << phi_x_error.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_y error range: [" << phi_y_error.min(0) << ", "
                   << phi_y_error.max(0) << "]" << std::endl;
}

void applyNestedL1InitializationScratchRepair(
    const NestedInitValidationContext& validation_ctx,
    NestedL1InitializationValidationScratch& scratch)
{
    if (validation_ctx.c3_use_interface == 1) {
        amrex::Print() << "\n    Using interface-based C3-repair (5-point reconstruction)" << std::endl;
        validation_ctx.apply_nested_c3_repair_pair(
            scratch.eta_prepared_temp, scratch.phi_prepared_temp, 1);
    } else {
        amrex::Print() << "    Using cell-center C3-repair (original)" << std::endl;
        crescent::nested::c3::applyC3SeamRepair(
            scratch.eta_prepared_temp,
            validation_ctx.geom,
            validation_ctx.grids,
            validation_ctx.domain,
            validation_ctx.c3_buffer_width,
            1);
        crescent::nested::c3::applyC3SeamRepair(
            scratch.phi_prepared_temp,
            validation_ctx.geom,
            validation_ctx.grids,
            validation_ctx.domain,
            validation_ctx.c3_buffer_width,
            1);
    }
}

void runNestedL1BoundaryRepairDiagnostic(
    const NestedInitValidationContext& validation_ctx,
    NestedL1InitializationValidationScratch& scratch)
{
    const auto before = collectNestedL1BoundaryDiagnostic(
        validation_ctx.grids, validation_ctx.geom, scratch.phi_prepared_temp);
    printNestedL1BoundaryDiagnosticBefore(before);
    applyNestedL1InitializationScratchRepair(validation_ctx, scratch);
    const auto after = collectNestedL1BoundaryDiagnostic(
        validation_ctx.grids, validation_ctx.geom, scratch.phi_prepared_temp);
    printNestedL1BoundaryDiagnosticAfter(before, after);
}

void computeNestedL1PreparedScratchDerivatives(
    const NestedInitValidationContext& validation_ctx,
    NestedL1InitializationValidationScratch& scratch)
{
    scratch.eta_prepared_x_temp.define(validation_ctx.grids, validation_ctx.dmap, 1, 0);
    scratch.eta_prepared_y_temp.define(validation_ctx.grids, validation_ctx.dmap, 1, 0);
    scratch.phi_prepared_x_temp.define(validation_ctx.grids, validation_ctx.dmap, 1, 0);
    scratch.phi_prepared_y_temp.define(validation_ctx.grids, validation_ctx.dmap, 1, 0);

    // Init-time local scratch for directional derivative spectra; hot-path
    // callers use the persistent SpectralWorkspace buffers instead.
    auto const& [cba, cdm] =
        validation_ctx.level_wave_operator.fft->getSpectralDataLayout();
    amrex::cMultiFab input_spectral_scratch(cba, cdm, 1, 0);
    amrex::cMultiFab x_spectral_scratch(cba, cdm, 1, 0);
    amrex::cMultiFab y_spectral_scratch(cba, cdm, 1, 0);

    crescent::physics::computeSurfaceDerivatives(
        scratch.eta_prepared_x_temp, scratch.eta_prepared_y_temp,
        scratch.phi_prepared_x_temp, scratch.phi_prepared_y_temp,
        scratch.eta_prepared_temp, scratch.phi_prepared_temp,
        validation_ctx.level_wave_operator.fft, validation_ctx.geom.data(),
        validation_ctx.spectral_settings.wavenumbers.dk_x,
        validation_ctx.spectral_settings.wavenumbers.dk_y,
        input_spectral_scratch,
        x_spectral_scratch,
        y_spectral_scratch);
}

void logNestedL1DerivativePipelineDiagnostic(
    const NestedInitValidationContext& validation_ctx,
    const NestedL1InitializationValidationScratch& scratch)
{
    amrex::Print() << "\n=== CRESCENT Derivative Pipeline Diagnostic ===" << std::endl;
    amrex::Print() << "  METHOD: C3-Periodic Seam Repair" << std::endl;
    amrex::Print() << "  phi_interp (from L0):        [" << scratch.phi_core_interp_reference.min(0)
                   << ", " << scratch.phi_core_interp_reference.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_periodic (C3-repair):    [" << scratch.phi_prepared_temp.min(0)
                   << ", " << scratch.phi_prepared_temp.max(0) << "]" << std::endl;

    amrex::Print() << "  phi_x_interp (from L0):      [" << scratch.phi_x_core_interp.min(0)
                   << ", " << scratch.phi_x_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_periodic_x (FFT):        [" << scratch.phi_prepared_x_temp.min(0)
                   << ", " << scratch.phi_prepared_x_temp.max(0) << "]" << std::endl;

    amrex::Print() << "  phi_x_blended (final):       [" << validation_ctx.phi_x.min(0)
                   << ", " << validation_ctx.phi_x.max(0) << "]" << std::endl;
    amrex::Print() << "  ---" << std::endl;
    amrex::Print() << "  phi_y_interp (should be 0):  [" << scratch.phi_y_core_interp.min(0)
                   << ", " << scratch.phi_y_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_periodic_y (should be 0):[" << scratch.phi_prepared_y_temp.min(0)
                   << ", " << scratch.phi_prepared_y_temp.max(0) << "]" << std::endl;

    amrex::Print() << "  phi_y_blended (should be 0): [" << validation_ctx.phi_y.min(0)
                   << ", " << validation_ctx.phi_y.max(0) << "]" << std::endl;
    amrex::Print() << "===============================================\n" << std::endl;
}

void blendNestedL1PreparedDiagnosticDerivatives(
    const NestedInitValidationContext& validation_ctx,
    NestedL1InitializationValidationScratch& scratch)
{
    scratch.phi_x_blended_temp.define(validation_ctx.grids, validation_ctx.dmap, 1, 0);
    scratch.phi_y_blended_temp.define(validation_ctx.grids, validation_ctx.dmap, 1, 0);

    for (MFIter mfi(scratch.phi_x_blended_temp); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& phi_x_out = scratch.phi_x_blended_temp.array(mfi);
        auto const& phi_x_per = scratch.phi_prepared_x_temp.const_array(mfi);
        auto const& phi_x_core_interp = scratch.phi_x_core_interp.const_array(mfi);
        auto const& sig = validation_ctx.sigma_blend.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real one_minus_sigma = 1.0 - sig(i, j, k);
            phi_x_out(i, j, k) = one_minus_sigma * phi_x_per(i, j, k)
                               + sig(i, j, k) * phi_x_core_interp(i, j, k);
        });
    }

    for (MFIter mfi(scratch.phi_y_blended_temp); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& phi_y_out = scratch.phi_y_blended_temp.array(mfi);
        auto const& phi_y_per = scratch.phi_prepared_y_temp.const_array(mfi);
        auto const& phi_y_core_interp = scratch.phi_y_core_interp.const_array(mfi);
        auto const& sig = validation_ctx.sigma_blend.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real one_minus_sigma = 1.0 - sig(i, j, k);
            phi_y_out(i, j, k) = one_minus_sigma * phi_y_per(i, j, k)
                               + sig(i, j, k) * phi_y_core_interp(i, j, k);
        });
    }

    amrex::Print() << "  Blended derivatives computed:" << std::endl;
    amrex::Print() << "    phi_x_blended: [" << scratch.phi_x_blended_temp.min(0)
                   << ", " << scratch.phi_x_blended_temp.max(0) << "]" << std::endl;
    amrex::Print() << "    phi_y_blended: [" << scratch.phi_y_blended_temp.min(0)
                   << ", " << scratch.phi_y_blended_temp.max(0) << "]" << std::endl;
}

void prepareNestedL1ValidationPeriodicCopies(
    const NestedInitValidationContext& validation_ctx,
    MultiFab& eta_before,
    MultiFab& phi_before,
    MultiFab& eta_periodic,
    MultiFab& phi_periodic)
{
    MultiFab::Copy(eta_before, validation_ctx.eta, 0, 0, 1, 0);
    MultiFab::Copy(phi_before, validation_ctx.phi, 0, 0, 1, 0);

    amrex::Print() << "\nApplying C3-repair to L1 domain (makes L1 boundaries periodic)..." << std::endl;
    amrex::Print() << "  L1 BoxArray: " << validation_ctx.grids << std::endl;
    amrex::Print() << "  Before C3-repair: eta range [" << validation_ctx.eta.min(0)
                   << ", " << validation_ctx.eta.max(0) << "]" << std::endl;
    amrex::Print() << "  Before C3-repair: phi range [" << validation_ctx.phi.min(0)
                   << ", " << validation_ctx.phi.max(0) << "]" << std::endl;

    MultiFab::Copy(eta_periodic, validation_ctx.eta, 0, 0, 1, 0);
    MultiFab::Copy(phi_periodic, validation_ctx.phi, 0, 0, 1, 0);

    if (validation_ctx.c3_use_interface == 1) {
        validation_ctx.apply_nested_c3_repair_pair(eta_periodic, phi_periodic, 0);
    } else {
        crescent::nested::c3::applyC3SeamRepair(
            eta_periodic, validation_ctx.geom, validation_ctx.grids,
            validation_ctx.c3_buffer_width, 0);
        crescent::nested::c3::applyC3SeamRepair(
            phi_periodic, validation_ctx.geom, validation_ctx.grids,
            validation_ctx.c3_buffer_width, 0);
    }

    amrex::Print() << "  L1 grid with periodic boundaries (ready for FFT):" << std::endl;
    amrex::Print() << "    eta range [" << eta_periodic.min(0) << ", "
                   << eta_periodic.max(0) << "]" << std::endl;
    amrex::Print() << "    phi range [" << phi_periodic.min(0) << ", "
                   << phi_periodic.max(0) << "]" << std::endl;
}

void computeAndLogNestedL1ValidationDerivatives(
    const NestedInitValidationContext& validation_ctx,
    const NestedL1InitializationValidationScratch& scratch,
    MultiFab& eta_x_CRESCENT,
    MultiFab& eta_y_CRESCENT,
    MultiFab& phi_x_CRESCENT,
    MultiFab& phi_y_CRESCENT,
    MultiFab& phi_z_CRESCENT)
{
    amrex::Print() << "\n=== Testing NEW L1-based derivative functions ===" << std::endl;
    crescent::nested::derivative_execution::computePreparedPeriodicDerivativeBundle(
        {
            validation_ctx.level_wave_operator,
            validation_ctx.geom,
            {
                scratch.eta_prepared_temp,
                scratch.phi_prepared_temp
            },
            {
                scratch.eta_prepared_temp,
                scratch.phi_prepared_temp
            },
            validation_ctx.box_S,
            validation_ctx.linear_derivative_spectral_filter_mode,
            validation_ctx.verbosity,
            validation_ctx.verbosity,
            {
                eta_x_CRESCENT,
                eta_y_CRESCENT,
                phi_x_CRESCENT,
                phi_y_CRESCENT,
                phi_z_CRESCENT
            },
            validation_ctx.wavenumber_pow,
            validation_ctx.eta_pow,
            validation_ctx.phi_expansion
        });

    amrex::Print() << "CRESCENT derivatives computed:" << std::endl;
    amrex::Print() << "  eta_x range: [" << eta_x_CRESCENT.min(0) << ", "
                   << eta_x_CRESCENT.max(0) << "]" << std::endl;
    amrex::Print() << "  eta_y range: [" << eta_y_CRESCENT.min(0) << ", "
                   << eta_y_CRESCENT.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_x range: [" << phi_x_CRESCENT.min(0) << ", "
                   << phi_x_CRESCENT.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_y range: [" << phi_y_CRESCENT.min(0) << ", "
                   << phi_y_CRESCENT.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_z range: [" << phi_z_CRESCENT.min(0) << ", "
                   << phi_z_CRESCENT.max(0) << "]" << std::endl;
}

void writeNestedL1PhiZDiagnostic(
    const NestedInitValidationContext& validation_ctx,
    NestedL1InitializationValidationScratch& scratch,
    const MultiFab& eta_x_CRESCENT,
    const MultiFab& eta_y_CRESCENT,
    const MultiFab& phi_x_CRESCENT,
    const MultiFab& phi_y_CRESCENT)
{
    amrex::Print() << "\n=== Computing phi_z diagnostic ===" << std::endl;

    constexpr int ng_interp = 2;

    MultiFab eta_pow_diag(
        validation_ctx.grids, validation_ctx.dmap,
        validation_ctx.perturbation_order - 1, 0);
    MultiFab phi_expansion_diag(
        validation_ctx.grids, validation_ctx.dmap,
        validation_ctx.perturbation_order, 0);

    amrex::Print() << "  phi[lev] (original): [" << validation_ctx.phi.min(0) << ", "
                   << validation_ctx.phi.max(0) << "]" << std::endl;
    amrex::Print() << "  phi_prepared_temp (C3-repaired): [" << scratch.phi_prepared_temp.min(0)
                   << ", " << scratch.phi_prepared_temp.max(0) << "]" << std::endl;

    MultiFab phi_z_nested_diag(validation_ctx.grids, validation_ctx.dmap, 1, 0);
    crescent::dno::computeNestedL1DNO(
        crescent::dno::makeSurfaceStateDNOContext(
            validation_ctx.level_wave_operator,
            validation_ctx.verbosity,
            {
                scratch.eta_prepared_temp,
                scratch.phi_prepared_temp
            },
            validation_ctx.wavenumber_pow,
            eta_pow_diag,
            phi_expansion_diag,
            phi_z_nested_diag));

    amrex::Print() << "  phi_z_L1 (from original eta[1]/phi[1]): ["
                   << phi_z_nested_diag.min(0) << ", " << phi_z_nested_diag.max(0) << "]" << std::endl;

    MultiFab phi_z_core_interp_l1(validation_ctx.grids, validation_ctx.dmap, 1, ng_interp);
    MultiFab phi_z_dummy(validation_ctx.grids, validation_ctx.dmap, 1, ng_interp);
    crescent::interpolation::interpolateParentToChild(
        validation_ctx.phi_z_core, validation_ctx.phi_z_core,
        phi_z_core_interp_l1, phi_z_dummy,
        validation_ctx.core_geom, validation_ctx.geom,
        validation_ctx.refine_ratio,
        scratch.coarse_interp_patch,
        validation_ctx.derivative_interpolation_order);
    amrex::Print() << "  phi_z_L0_interp: [" << phi_z_core_interp_l1.min(0) << ", "
                   << phi_z_core_interp_l1.max(0) << "]" << std::endl;

    MultiFab phi_z_blended(validation_ctx.grids, validation_ctx.dmap, 1, 0);
    for (MFIter mfi(phi_z_blended); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto phi_z_nested_arr = phi_z_nested_diag.const_array(mfi);
        auto phi_z_core_arr = phi_z_core_interp_l1.const_array(mfi);
        auto sigma_arr = validation_ctx.sigma_blend.const_array(mfi);
        auto phi_z_blend_arr = phi_z_blended.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Real sig = sigma_arr(i, j, k);
            phi_z_blend_arr(i, j, k) = (1.0 - sig) * phi_z_nested_arr(i, j, k)
                                     + sig * phi_z_core_arr(i, j, k);
        });
    }
    amrex::Print() << "  phi_z_blended: [" << phi_z_blended.min(0) << ", "
                   << phi_z_blended.max(0) << "]" << std::endl;

    crescent::nested::diagnostic::writeDiagnosticTecplot(
        validation_ctx.diagnostic_output_file,
        validation_ctx.geom,
        scratch.eta_core_interp, scratch.eta_prepared_temp,
        scratch.eta_x_core_interp, eta_x_CRESCENT,
        scratch.eta_y_core_interp, eta_y_CRESCENT,
        scratch.phi_core_interp, scratch.phi_prepared_temp,
        scratch.phi_x_core_interp, scratch.phi_prepared_x_temp, scratch.phi_x_blended_temp,
        scratch.phi_y_core_interp, scratch.phi_prepared_y_temp, scratch.phi_y_blended_temp,
        phi_z_nested_diag, phi_z_core_interp_l1, phi_z_blended,
        validation_ctx.root_ba, validation_ctx.root_dmap);
    amrex::Print() << "=== phi_z diagnostic complete ===" << std::endl;
}

} // namespace

namespace crescent::nested::validation {

void runInitializationValidation(
    const InitializationValidationContext& context)
{
    const int lev = context.lev;
    auto& solution = context.core.solution;
    const auto& runtime = context.core.runtime;
    const auto& wave = context.core.wave;
    const auto& nested = context.nested.config;
    AMREX_ALWAYS_ASSERT(nested.use_c3_repair == 1);
    if (runtime.verbosity < 2) {
        return;
    }

    const auto core_mesh = context.coreLevelMesh();
    const auto level_mesh = context.levelMesh();
    const auto refine_ratio = context.refineRatio();
    const auto& level_wave_operator = context.level_wave_operator;
    amrex::Print() << "\n=== CRESCENT_Derivatives Validation ===" << std::endl;

    NestedL1InitializationValidationScratch scratch;
    const BoxArray& ba_nested = solution.eta[lev].boxArray();
    const DistributionMapping& dm_nested = solution.eta[lev].DistributionMap();
    constexpr int ng_interp = 2;

    scratch.eta_core_interp.define(ba_nested, dm_nested, 1, ng_interp);
    scratch.phi_core_interp.define(ba_nested, dm_nested, 1, ng_interp);
    scratch.eta_x_core_interp.define(ba_nested, dm_nested, 1, ng_interp);
    scratch.eta_y_core_interp.define(ba_nested, dm_nested, 1, ng_interp);
    scratch.phi_x_core_interp.define(ba_nested, dm_nested, 1, ng_interp);
    scratch.phi_y_core_interp.define(ba_nested, dm_nested, 1, ng_interp);
    scratch.coarse_interp_patch.define(
        amrex::coarsen(ba_nested, refine_ratio),
        dm_nested,
        1,
        ng_interp);

    crescent::interpolation::interpolateParentToChild(
        solution.eta[0], solution.phi[0], scratch.eta_core_interp, scratch.phi_core_interp,
        core_mesh.geom, level_mesh.geom, refine_ratio,
        scratch.coarse_interp_patch, nested.state_interpolation_order);

    crescent::interpolation::interpolateParentToChild(
        solution.eta_x[0], solution.phi_x[0], scratch.eta_x_core_interp, scratch.phi_x_core_interp,
        core_mesh.geom, level_mesh.geom, refine_ratio,
        scratch.coarse_interp_patch, nested.derivative_interpolation_order);
    crescent::interpolation::interpolateParentToChild(
        solution.eta_y[0], solution.phi_y[0], scratch.eta_y_core_interp, scratch.phi_y_core_interp,
        core_mesh.geom, level_mesh.geom, refine_ratio,
        scratch.coarse_interp_patch, nested.derivative_interpolation_order);

    amrex::Print() << "    Interpolated L0 fields and derivatives:" << std::endl;
    amrex::Print() << "      eta_L0_interp range: [" << scratch.eta_core_interp.min(0)
                   << ", " << scratch.eta_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "      phi_L0_interp range: [" << scratch.phi_core_interp.min(0)
                   << ", " << scratch.phi_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "      eta_x_L0_interp range: [" << scratch.eta_x_core_interp.min(0)
                   << ", " << scratch.eta_x_core_interp.max(0) << "]" << std::endl;
    amrex::Print() << "      phi_x_L0_interp range: [" << scratch.phi_x_core_interp.min(0)
                   << ", " << scratch.phi_x_core_interp.max(0) << "]" << std::endl;

    amrex::Print() << "    Preparing initialization-only repaired copies for diagnostics" << std::endl;

    scratch.phi_core_interp_reference.define(ba_nested, dm_nested, 1, ng_interp);
    crescent::interpolation::interpolateParentToChild(
        solution.phi[0], solution.phi[0], scratch.phi_core_interp_reference, scratch.phi_core_interp_reference,
        core_mesh.geom, level_mesh.geom, refine_ratio,
        scratch.coarse_interp_patch, nested.state_interpolation_order);

    scratch.eta_prepared_temp.define(ba_nested, dm_nested, 1, 0);
    scratch.phi_prepared_temp.define(ba_nested, dm_nested, 1, 0);
    MultiFab::Copy(scratch.eta_prepared_temp, solution.eta[lev], 0, 0, 1, 0);
    MultiFab::Copy(scratch.phi_prepared_temp, solution.phi[lev], 0, 0, 1, 0);
    const auto spectral_settings = level_wave_operator.spectral_settings;

    const NestedInitValidationContext validation_ctx{
        wave.spectral.perturbation_order,
        runtime.verbosity,
        nested.c3_use_interface,
        nested.linear_derivative_spectral_filter_mode,
        nested.derivative_interpolation_order,
        spectral_settings,
        refine_ratio,
        level_mesh.grids,
        level_mesh.dmap,
        level_mesh.geom,
        core_mesh.geom,
        context.domain.box_S,
        context.root_grids,
        context.root_dmap,
        context.domain,
        context.c3_buffer_width,
        level_wave_operator,
        solution.eta[lev],
        solution.phi[lev],
        solution.phi_x[lev],
        solution.phi_y[lev],
        solution.wavenumber_pow[lev],
        solution.eta_pow[lev],
        solution.phi_expansion[lev],
        solution.phi_z[0],
        context.sigma_blend,
        context.diagnostic_output_file,
        [&context](amrex::MultiFab& eta_mf,
                   amrex::MultiFab& phi_mf,
                   int repair_verbosity) {
            crescent::nested::c3::applyRepairPair(
                context.lev, eta_mf, phi_mf,
                context.mesh.geom, context.mesh.grids,
                context.nested.config, context.nested.runtime.coupling_zones,
                repair_verbosity,
                context.nested.runtime.c3_repair_workspace[context.lev]);
        }
    };

    runNestedL1BoundaryRepairDiagnostic(validation_ctx, scratch);

    computeNestedL1PreparedScratchDerivatives(validation_ctx, scratch);
    logNestedL1DerivativePipelineDiagnostic(validation_ctx, scratch);
    blendNestedL1PreparedDiagnosticDerivatives(validation_ctx, scratch);

    MultiFab eta_x_CRESCENT(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab eta_y_CRESCENT(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab phi_x_CRESCENT(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab phi_y_CRESCENT(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab phi_z_CRESCENT(level_mesh.grids, level_mesh.dmap, 1, 0);

    MultiFab eta_before(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab phi_before(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab eta_periodic(level_mesh.grids, level_mesh.dmap, 1, 0);
    MultiFab phi_periodic(level_mesh.grids, level_mesh.dmap, 1, 0);

    prepareNestedL1ValidationPeriodicCopies(
        validation_ctx, eta_before, phi_before, eta_periodic, phi_periodic);
    writeNestedL1PeriodicitySlice(
        level_mesh.grids, level_mesh.geom,
        solution.eta[lev], solution.phi[lev], eta_before, phi_before);
    computeAndLogNestedL1ValidationDerivatives(
        validation_ctx, scratch,
        eta_x_CRESCENT, eta_y_CRESCENT, phi_x_CRESCENT, phi_y_CRESCENT, phi_z_CRESCENT);
    writeNestedL1PhiZDiagnostic(
        validation_ctx, scratch,
        eta_x_CRESCENT, eta_y_CRESCENT, phi_x_CRESCENT, phi_y_CRESCENT);
    logNestedL1InterpolatedDerivativeErrors(
        level_mesh.grids, level_mesh.dmap, scratch,
        eta_x_CRESCENT, eta_y_CRESCENT, phi_x_CRESCENT, phi_y_CRESCENT);

    amrex::Print() << "=== CRESCENT_Derivatives Validation Complete ===" << std::endl;
}

} // namespace crescent::nested::validation

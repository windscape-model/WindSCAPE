#include "CRESCENT_NestedCorrectionField.H"

#include <cmath>
#include <memory>

#include <AMReX_MultiFab.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_Dynamics.H"
#include "CRESCENT_Interpolation.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_NestedRK4.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_SpectralFilterExecution.H"
#include "CRESCENT_SpectralFilterKernels.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_SolverContext.H"
#include "shared/CRESCENT_LevelWaveOperatorView.H"

namespace crescent::nested::correction_field {

namespace {

using amrex::Box;
using amrex::GpuComplex;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;
using crescent::shared::LevelContext;

// In-place spectral multiply by sign |k|^power (Riesz symbol) on a parent-
// or child-grid MultiFab, optionally with the per-axis square two-thirds
// cut of the parent (spectral::dealiasing and the DNO recursion).
void applySpectralRieszMultiplier(
    MultiFab& mf,
    std::unique_ptr<amrex::FFT::R2C<Real>>& fft,
    const amrex::GeometryData& geomdata,
    Real sign,
    int power,
    bool two_thirds_cutoff)
{
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(mf, geomdata);
    const int nx = spectral_grid.nx;
    const int ny = spectral_grid.ny;
    const Real dk_x = spectral_grid.dk_x;
    const Real dk_y = spectral_grid.dk_y;
    const Real norm = Real(spectral_grid.nx) * Real(spectral_grid.ny);
    const int cutoff_x =
        static_cast<int>(nx * crescent::constants::ONE_THIRD_CUTOFF);
    const int cutoff_y =
        static_cast<int>(ny * crescent::constants::ONE_THIRD_CUTOFF);

    const auto multiplier =
        [=] AMREX_GPU_DEVICE(int i, int j, int, GpuComplex<Real>& fft_c) noexcept {
            const Real k_mag =
                crescent::spectral::filter_kernels::spectralMagnitude(
                    i, j, ny, dk_x, dk_y);
            Real weight = sign;
            for (int p = 0; p < power; ++p) {
                weight *= k_mag;
            }
            if (two_thirds_cutoff &&
                crescent::spectral::filter_kernels::spectralModeExceedsCutoff(
                    crescent::spectral::filter_kernels::spectralKxIndex(i, nx),
                    crescent::spectral::filter_kernels::spectralKyIndex(j, ny),
                    cutoff_x, cutoff_y)) {
                weight = 0.0;
            }
            fft_c = fft_c * (weight / norm);
        };

    crescent::spectral::applyInPlaceSpectralFilter(mf, fft, multiplier);
}

// out = a * b (component 0 each), valid region.
void multiplyInto(MultiFab& out, const MultiFab& a, const MultiFab& b)
{
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto o = out.array(mfi);
        auto aa = a.const_array(mfi);
        auto bb = b.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            o(i, j, k) = aa(i, j, k) * bb(i, j, k);
        });
    }
}

struct CorrectionStageContext
{
    int lev;
    crescent::shared::LevelContext& level_context;
    crescent::nested::CorrectionFieldScratch& cf;
    const MultiFab& chi_source_sigma; // sigma; chi = 1 - sigma
    crescent::shared::LevelWaveOperatorView child_operator;
};

CorrectionStageContext makeCorrectionStageContext(
    crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    return {
        lev,
        level_context,
        level_context.nestedRuntime().correction_field,
        level_context.nestedRuntime().coupling_zones.sigma_blend[lev],
        level_context.waveOperator()
    };
}

// out = chi * (total - baseline)
void formTaperedCorrection(
    MultiFab& out,
    const MultiFab& total,
    const MultiFab& baseline,
    const MultiFab& sigma)
{
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto o = out.array(mfi);
        auto t = total.const_array(mfi);
        auto b = baseline.const_array(mfi);
        auto s = sigma.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            o(i, j, k) = (1.0 - s(i, j, k)) * (t(i, j, k) - b(i, j, k));
        });
    }
}

void interpolateToChild(
    const crescent::shared::LevelContext& level_context,
    const MultiFab& coarse,
    MultiFab& fine)
{
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;
    crescent::interpolation::interpolateScalarParentToChild(
        coarse,
        fine,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        coupling_derivatives.core_interp_patch_scratch,
        level_context.nestedConfig().state_interpolation_order);
}

// baseline_on_child = I( saved_L0 + fac * rhs_L0[comp] ); fac==0 -> I(saved).
void reconstructBaselineOnChild(
    const crescent::shared::LevelContext& level_context,
    const MultiFab& saved_core,
    const MultiFab& rhs_core,
    int rhs_comp,
    Real fac,
    MultiFab& core_scratch,
    MultiFab& baseline_on_child)
{
    for (MFIter mfi(core_scratch, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto o = core_scratch.array(mfi);
        auto s = saved_core.const_array(mfi);
        auto r = rhs_core.const_array(mfi);
        const int comp = rhs_comp;
        const Real f = fac;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            o(i, j, k) = s(i, j, k) + (f != 0.0 ? f * r(i, j, k, comp) : 0.0);
        });
    }
    interpolateToChild(level_context, core_scratch, baseline_on_child);
}

// Order-2 bg-delta cross terms of dw (spec §4). Requires M >= 2.
// Uses cf.cross_a / cf.cross_b as working buffers; adds into cf.dw.
void accumulateCrossTerms(const CorrectionStageContext& ctx)
{
    auto& cf = ctx.cf;
    auto& fft = ctx.child_operator.fft;
    const auto& geomdata = ctx.child_operator.geomdata;

    // B1 = -eta0 * lap(T(dp)): the same coefficient cutoff T on phi^(0) that
    // the DNO recurrence applies before vertical-velocity accumulation.
    MultiFab::Copy(cf.cross_a, cf.dp, 0, 0, 1, 0);
    applySpectralRieszMultiplier(cf.cross_a, fft, geomdata, -1.0, 2, true);
    for (MFIter mfi(cf.dw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto dw = cf.dw.array(mfi);
        auto eta0 = cf.bg_on_l1.eta.const_array(mfi);
        auto lap_dp = cf.cross_a.const_array(mfi);
        auto lap0 = cf.bg_on_l1.lap_phi.const_array(mfi);
        auto de = cf.de.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // B1 + B2 = -(eta0 * lap(T(dp)) + de * I(lap(T0(phi0))))
            dw(i, j, k) -= eta0(i, j, k) * lap_dp(i, j, k)
                         + de(i, j, k) * lap0(i, j, k);
        });
    }

    // A1 = -IFFT(|k| Cut23[FFT( chi * eta0 * IFFT(|k| dp_hat) )]):
    // cross_a = IFFT(|k| dp_hat); then taper the bg product; then outer op.
    MultiFab::Copy(cf.cross_a, cf.dp, 0, 0, 1, 0);
    applySpectralRieszMultiplier(cf.cross_a, fft, geomdata, 1.0, 1, false);
    for (MFIter mfi(cf.cross_a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto a = cf.cross_a.array(mfi);
        auto eta0 = cf.bg_on_l1.eta.const_array(mfi);
        auto s = ctx.chi_source_sigma.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            a(i, j, k) *= (1.0 - s(i, j, k)) * eta0(i, j, k);
        });
    }
    applySpectralRieszMultiplier(cf.cross_a, fft, geomdata, -1.0, 1, true);

    // A2 = -IFFT(|k| Cut23[FFT( de * I(w_lin0) )])
    multiplyInto(cf.cross_b, cf.de, cf.bg_on_l1.w_lin);
    applySpectralRieszMultiplier(cf.cross_b, fft, geomdata, -1.0, 1, true);

    MultiFab::Add(cf.dw, cf.cross_a, 0, 0, 1, 0);
    MultiFab::Add(cf.dw, cf.cross_b, 0, 0, 1, 0);
}

void interpolateStageBackground(
    const CorrectionStageContext& ctx,
    const crescent::nested::CoreStageSnapshotFields& capture)
{
    auto& bg = ctx.cf.bg_on_l1;
    interpolateToChild(ctx.level_context, capture.eta, bg.eta);
    interpolateToChild(ctx.level_context, capture.eta_x, bg.eta_x);
    interpolateToChild(ctx.level_context, capture.eta_y, bg.eta_y);
    interpolateToChild(ctx.level_context, capture.phi_x, bg.phi_x);
    interpolateToChild(ctx.level_context, capture.phi_y, bg.phi_y);
    interpolateToChild(ctx.level_context, capture.phi_z, bg.phi_z);
    interpolateToChild(ctx.level_context, capture.w_lin, bg.w_lin);
    interpolateToChild(ctx.level_context, capture.lap_phi, bg.lap_phi);
    interpolateToChild(ctx.level_context, capture.m0_d, bg.m0_d);
    interpolateToChild(ctx.level_context, capture.w0sq_d, bg.w0sq_d);
}

void dealiasOnChild(const CorrectionStageContext& ctx, MultiFab& mf)
{
    auto& workspace = ctx.child_operator.spectral_workspace;
    crescent::spectral::dealiasing(
        mf,
        ctx.child_operator.fft,
        ctx.child_operator.geomdata,
        workspace.dealiasing_scratch[ctx.lev]);
}

void applyCorrectionTaper(
    const CorrectionStageContext& ctx,
    MultiFab& field)
{
    for (MFIter mfi(field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto values = field.array(mfi);
        auto sigma = ctx.chi_source_sigma.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            values(i, j, k) *= (1.0 - sigma(i, j, k));
        });
    }
}

// Parent-parity final two-thirds projection of a correction DNO field. The
// background-cross terms are not periodic until they are compactly
// supported: taper before the FFT projection, then taper again because the
// projection is nonlocal and can repopulate the exterior buffer.
void finalizeCompactDNOField(
    const CorrectionStageContext& ctx,
    MultiFab& field)
{
    applyCorrectionTaper(ctx, field);
    dealiasOnChild(ctx, field);
    applyCorrectionTaper(ctx, field);
}

// Added-band projector: zero every child mode
// inside the parent's retained square band, |kx| <= k_d and |ky| <= k_d with
// k_d = int(nx0/3) * dk_x0, so that only the modes the parent cannot represent
// remain. The input must be compactly supported or a child Fourier series
// (periodicity discipline); the output is a child Fourier series.
void applyAddedBandProjector(
    const CorrectionStageContext& ctx,
    MultiFab& field)
{
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(field, ctx.child_operator.geomdata);
    const int nx = spectral_grid.nx;
    const int ny = spectral_grid.ny;
    const Real norm = Real(nx) * Real(ny);

    const amrex::Geometry& core_geom = ctx.level_context.coreLevelGeometry();
    const int nx0 = core_geom.Domain().length(0);
    const int ny0 = core_geom.Domain().length(1);
    const Real dk_x0 = 2.0 * crescent::constants::PI / core_geom.ProbLength(0);
    const Real dk_y0 = 2.0 * crescent::constants::PI / core_geom.ProbLength(1);
    const Real k_shared_x =
        static_cast<int>(nx0 * crescent::constants::ONE_THIRD_CUTOFF) * dk_x0;
    const Real k_shared_y =
        static_cast<int>(ny0 * crescent::constants::ONE_THIRD_CUTOFF) * dk_y0;
    const int cutoff_x =
        static_cast<int>(std::floor(k_shared_x / spectral_grid.dk_x + 1.0e-9));
    const int cutoff_y =
        static_cast<int>(std::floor(k_shared_y / spectral_grid.dk_y + 1.0e-9));

    const auto projector =
        [=] AMREX_GPU_DEVICE(int i, int j, int, GpuComplex<Real>& fft_c) noexcept {
            const int kx =
                crescent::spectral::filter_kernels::spectralKxIndex(i, nx);
            const int ky =
                crescent::spectral::filter_kernels::spectralKyIndex(j, ny);
            const int akx = (kx < 0) ? -kx : kx;
            const int aky = (ky < 0) ? -ky : ky;
            if (akx <= cutoff_x && aky <= cutoff_y) {
                fft_c = GpuComplex<Real>(0.0, 0.0);
            } else {
                fft_c = fft_c / norm;
            }
        };

    crescent::spectral::applyInPlaceSpectralFilter(field, ctx.child_operator.fft, projector);
}

} // namespace

namespace {

// Fill a parent-grid background snapshot from a parent state pair and its
// derivative bundle: copies, the two parent-grid spectral quantities
// w_lin = IFFT(|k| phi_hat) and lap_phi = IFFT(-|k|^2 phi_hat), and the
// parent-dealiased products m0_d, w0sq_d (spec §4).
void fillCoreBackgroundSnapshot(
    const crescent::shared::LevelWaveOperatorView& core_operator,
    const MultiFab& eta0,
    const MultiFab& phi0,
    const MultiFab& eta_x0,
    const MultiFab& eta_y0,
    const MultiFab& phi_x0,
    const MultiFab& phi_y0,
    const MultiFab& phi_z0,
    crescent::nested::CoreStageSnapshotFields& snap)
{
    MultiFab::Copy(snap.eta, eta0, 0, 0, 1, 0);
    MultiFab::Copy(snap.phi, phi0, 0, 0, 1, 0);
    MultiFab::Copy(snap.eta_x, eta_x0, 0, 0, 1, 0);
    MultiFab::Copy(snap.eta_y, eta_y0, 0, 0, 1, 0);
    MultiFab::Copy(snap.phi_x, phi_x0, 0, 0, 1, 0);
    MultiFab::Copy(snap.phi_y, phi_y0, 0, 0, 1, 0);
    MultiFab::Copy(snap.phi_z, phi_z0, 0, 0, 1, 0);

    auto& fft0 = core_operator.fft;
    const auto& geomdata0 = core_operator.geomdata;

    // w_lin uses the unprojected surface potential because it appears inside
    // the phi^(1) recurrence. The W1 B term uses the two-thirds-projected
    // phi^(0), matching vertical-velocity accumulation.
    MultiFab::Copy(snap.w_lin, snap.phi, 0, 0, 1, 0);
    applySpectralRieszMultiplier(snap.w_lin, fft0, geomdata0, 1.0, 1, false);
    MultiFab::Copy(snap.lap_phi, snap.phi, 0, 0, 1, 0);
    applySpectralRieszMultiplier(snap.lap_phi, fft0, geomdata0, -1.0, 2, true);

    // Parent-dealiased background products.
    for (MFIter mfi(snap.m0_d, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto m0 = snap.m0_d.array(mfi);
        auto wsq = snap.w0sq_d.array(mfi);
        auto ex = snap.eta_x.const_array(mfi);
        auto ey = snap.eta_y.const_array(mfi);
        auto w = snap.phi_z.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            m0(i, j, k) = ex(i, j, k) * ex(i, j, k) + ey(i, j, k) * ey(i, j, k);
            wsq(i, j, k) = w(i, j, k) * w(i, j, k);
        });
    }
    auto& workspace = core_operator.spectral_workspace;
    crescent::spectral::dealiasing(
        snap.m0_d, fft0, geomdata0, workspace.dealiasing_scratch[0]);
    crescent::spectral::dealiasing(
        snap.w0sq_d, fft0, geomdata0, workspace.dealiasing_scratch[0]);
}

} // namespace

void captureCoreStage(
    const crescent::shared::LevelContext& core_level_context,
    int rk_stage)
{
    if (core_level_context.level() != 0) {
        return;
    }
    auto& solution = core_level_context.solutionState();
    auto& cf = core_level_context.nestedRuntime().correction_field;
    if (!cf.core_stage_defined) {
        return;
    }
    const auto core_operator = core_level_context.waveOperator();
    fillCoreBackgroundSnapshot(
        core_operator,
        solution.eta[0], solution.phi[0],
        solution.eta_x[0], solution.eta_y[0],
        solution.phi_x[0], solution.phi_y[0], solution.phi_z[0],
        cf.core_stage[rk_stage]);
}

namespace {

// The steps of one correction-field RK4 stage (runL1CorrectionFieldStage).

// Step 1: the pre-filter parent stage baseline u0 = I(saved + fac * rhs) on
// the child (spec section 1) and the tapered correction (de, dp) = chi (u - u0).
void formStageCorrection(const CorrectionStageContext& ctx, int irk, Real fac)
{
    auto& level_context = ctx.level_context;
    const int lev = ctx.lev;
    auto& solution = level_context.solutionState();
    auto& cf = ctx.cf;
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;

    reconstructBaselineOnChild(
        level_context,
        solution.eta_saved[0],
        solution.rhs_eta[0],
        irk > 0 ? irk - 1 : 0,
        irk > 0 ? fac : 0.0,
        coupling_derivatives.core_stage_derivative_scratch,
        cf.cross_a);
    formTaperedCorrection(
        cf.de, solution.eta[lev], cf.cross_a, ctx.chi_source_sigma);
    reconstructBaselineOnChild(
        level_context,
        solution.phi_saved[0],
        solution.rhs_phi[0],
        irk > 0 ? irk - 1 : 0,
        irk > 0 ? fac : 0.0,
        coupling_derivatives.core_stage_derivative_scratch,
        cf.cross_a);
    formTaperedCorrection(
        cf.dp, solution.phi[lev], cf.cross_a, ctx.chi_source_sigma);
}

// Steps 3-4: the spectral derivatives of the correction (into the solution
// bundle, which the correction-field coupling does not otherwise use) and its
// vertical velocity to full order (M = 1 gives the plain |k| chain).
void differentiateCorrection(const CorrectionStageContext& ctx)
{
    const int lev = ctx.lev;
    auto& solution = ctx.level_context.solutionState();
    auto& cf = ctx.cf;

    auto& workspace = ctx.child_operator.spectral_workspace;
    crescent::physics::computeSurfaceDerivatives(
        solution.eta_x[lev],
        solution.eta_y[lev],
        solution.phi_x[lev],
        solution.phi_y[lev],
        cf.de,
        cf.dp,
        ctx.child_operator.fft,
        ctx.child_operator.geomdata,
        ctx.child_operator.spectral_settings.wavenumbers.dk_x,
        ctx.child_operator.spectral_settings.wavenumbers.dk_y,
        workspace.directional_input_spectral[lev],
        workspace.directional_x_spectral[lev],
        workspace.directional_y_spectral[lev]);

    crescent::dno::computeNestedL1DNO(
        ctx.child_operator,
        cf.de,
        cf.dp,
        solution.wavenumber_pow[lev],
        solution.eta_pow[lev],
        solution.phi_expansion[lev],
        cf.dw,
        false);
}

// Step 5, order M = 1 (spec section 5): with the child-assembled background,
// F_eta = w0 + dw and F_psi = -g (eta0 + de) - p1.
void assembleLinearStageRHS(const CorrectionStageContext& ctx, int irk)
{
    const int lev = ctx.lev;
    auto& solution = ctx.level_context.solutionState();
    auto& cf = ctx.cf;

    finalizeCompactDNOField(ctx, cf.dw);

    amrex::MultiFab rhs_eta_stage(
        solution.rhs_eta[lev], amrex::make_alias, irk, 1);
    amrex::MultiFab rhs_phi_stage(
        solution.rhs_phi[lev], amrex::make_alias, irk, 1);
    for (MFIter mfi(rhs_eta_stage, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto re = rhs_eta_stage.array(mfi);
        auto rp = rhs_phi_stage.array(mfi);
        auto dw = cf.dw.const_array(mfi);
        auto de = cf.de.const_array(mfi);
        auto e0 = cf.bg_on_l1.eta.const_array(mfi);
        auto w0 = cf.bg_on_l1.phi_z.const_array(mfi);
        auto p1 = solution.pressure[lev].const_array(mfi);
        const Real g = crescent::constants::GRAVITY;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            re(i, j, k) = w0(i, j, k) + dw(i, j, k);
            rp(i, j, k) = -g * (e0(i, j, k) + de(i, j, k))
                          - p1(i, j, k);
        });
    }
}

// Step 7: the metric increment dm_raw = 2 grad(eta0).grad(de) + |grad(de)|^2
// and its tapered, child-dealiased copy dm. The background cross term is
// non-periodic after the spectral derivative of de fills the box, so the
// taper restores compact support before the transform.
void formMetricIncrement(const CorrectionStageContext& ctx)
{
    const int lev = ctx.lev;
    auto& solution = ctx.level_context.solutionState();
    auto& cf = ctx.cf;

    for (MFIter mfi(cf.dm_raw, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto dmr = cf.dm_raw.array(mfi);
        auto e0x = cf.bg_on_l1.eta_x.const_array(mfi);
        auto e0y = cf.bg_on_l1.eta_y.const_array(mfi);
        auto dex = solution.eta_x[lev].const_array(mfi);
        auto dey = solution.eta_y[lev].const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            dmr(i, j, k) = 2.0 * (e0x(i, j, k) * dex(i, j, k)
                                  + e0y(i, j, k) * dey(i, j, k))
                           + dex(i, j, k) * dex(i, j, k)
                           + dey(i, j, k) * dey(i, j, k);
        });
    }
    MultiFab::Copy(cf.dm, cf.dm_raw, 0, 0, 1, 0);
    applyCorrectionTaper(ctx, cf.dm);
    dealiasOnChild(ctx, cf.dm);
}

// Step 8: the dealiased w-square increment term_b = Dealias[2 w0 dw + dw^2].
void formWSquareIncrement(const CorrectionStageContext& ctx)
{
    auto& cf = ctx.cf;

    for (MFIter mfi(cf.term_b, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto tb = cf.term_b.array(mfi);
        auto w0 = cf.bg_on_l1.phi_z.const_array(mfi);
        auto dw = cf.dw.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tb(i, j, k) = 2.0 * w0(i, j, k) * dw(i, j, k)
                          + dw(i, j, k) * dw(i, j, k);
        });
    }
    dealiasOnChild(ctx, cf.term_b);
}

// Step 9: dm_w2 = Dealias[chi * dm_raw * (I(w0sq_d) + 2 w0 dw + dw^2)] in
// cross_b.
void formMetricWSquareTerm(const CorrectionStageContext& ctx)
{
    auto& cf = ctx.cf;

    for (MFIter mfi(cf.cross_b, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto cb = cf.cross_b.array(mfi);
        auto dmr = cf.dm_raw.const_array(mfi);
        auto w0sq = cf.bg_on_l1.w0sq_d.const_array(mfi);
        auto w0 = cf.bg_on_l1.phi_z.const_array(mfi);
        auto dw = cf.dw.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            cb(i, j, k) = dmr(i, j, k)
                          * (w0sq(i, j, k)
                             + 2.0 * w0(i, j, k) * dw(i, j, k)
                             + dw(i, j, k) * dw(i, j, k));
        });
    }
    applyCorrectionTaper(ctx, cf.cross_b);
    dealiasOnChild(ctx, cf.cross_b);
}

// Step 10: the tendency increment DeltaF (spec section 3) in term_a (eta)
// and term_b (psi).
void assembleTendencyIncrement(const CorrectionStageContext& ctx)
{
    const int lev = ctx.lev;
    auto& solution = ctx.level_context.solutionState();
    auto& cf = ctx.cf;

    for (MFIter mfi(cf.term_a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto ta = cf.term_a.array(mfi);
        auto tb = cf.term_b.array(mfi);
        auto de = cf.de.const_array(mfi);
        auto dex = solution.eta_x[lev].const_array(mfi);
        auto dey = solution.eta_y[lev].const_array(mfi);
        auto dpx = solution.phi_x[lev].const_array(mfi);
        auto dpy = solution.phi_y[lev].const_array(mfi);
        auto e0x = cf.bg_on_l1.eta_x.const_array(mfi);
        auto e0y = cf.bg_on_l1.eta_y.const_array(mfi);
        auto p0x = cf.bg_on_l1.phi_x.const_array(mfi);
        auto p0y = cf.bg_on_l1.phi_y.const_array(mfi);
        auto w0 = cf.bg_on_l1.phi_z.const_array(mfi);
        auto m0d = cf.bg_on_l1.m0_d.const_array(mfi);
        auto dm = cf.dm.const_array(mfi);
        auto dw = cf.dw.const_array(mfi);
        auto dm_w2 = cf.cross_b.const_array(mfi);
        const Real g = crescent::constants::GRAVITY;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            ta(i, j, k) =
                -(p0x(i, j, k) * dex(i, j, k) + p0y(i, j, k) * dey(i, j, k))
                - (dpx(i, j, k) * e0x(i, j, k) + dpy(i, j, k) * e0y(i, j, k))
                - (dpx(i, j, k) * dex(i, j, k) + dpy(i, j, k) * dey(i, j, k))
                + (1.0 + m0d(i, j, k) + dm(i, j, k)) * dw(i, j, k)
                + dm(i, j, k) * w0(i, j, k);
            // tb currently holds Dealias[2 w0 dwt + dwt^2]; overwrite with
            // the full DeltaF_psi using it plus the precomputed dm_w2.
            tb(i, j, k) =
                -g * de(i, j, k)
                - (p0x(i, j, k) * dpx(i, j, k) + p0y(i, j, k) * dpy(i, j, k))
                - 0.5 * (dpx(i, j, k) * dpx(i, j, k)
                         + dpy(i, j, k) * dpy(i, j, k))
                + 0.5 * (1.0 + m0d(i, j, k)) * tb(i, j, k)
                + 0.5 * dm_w2(i, j, k);
        });
    }
}

// Step 11: the final taper and dealias of DeltaF, at the points where the
// parent dealiases its tendency.
void taperAndDealiasIncrement(const CorrectionStageContext& ctx)
{
    auto& cf = ctx.cf;

    for (MFIter mfi(cf.term_a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto ta = cf.term_a.array(mfi);
        auto tb = cf.term_b.array(mfi);
        auto s = ctx.chi_source_sigma.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real chi = 1.0 - s(i, j, k);
            ta(i, j, k) *= chi;
            tb(i, j, k) *= chi;
        });
    }
    dealiasOnChild(ctx, cf.term_a);
    dealiasOnChild(ctx, cf.term_b);
}

struct StageAdjustment
{
    bool active = false;
    Real weight = 1.0;
};

// The Dommermuth adjustment of the stage: the weight W of the nonlinear
// remainder at the stage time and, while W < 1, the linear part of the
// compact correction, built through the same taper and projection chain as
// the corresponding terms of DeltaF, so that the stage uses
// L(delta) + W * (DeltaF - L(delta)) without ramping pressure. The linear
// parts go to adjustment_w_lin_scratch (eta) and cf.dm_raw (phi state).
StageAdjustment prepareStageAdjustment(
    const CorrectionStageContext& ctx,
    Real stage_time)
{
    auto& level_context = ctx.level_context;
    auto& cf = ctx.cf;
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;

    const auto& wave = level_context.core().wave;
    Real adjustment_weight = 1.0;
    bool adjustment_active = false;
    if (wave.spectral.adjustment_periods > 0.0) {
        const Real t_peak = 2.0 * crescent::constants::PI /
            wave.initialization.peak_frequency;
        adjustment_weight = crescent::dynamics::computeAdjustmentRampWeight(
            stage_time,
            wave.spectral.adjustment_periods * t_peak);
        adjustment_active = adjustment_weight < 1.0;
    }

    auto& linear_delta_eta = coupling_derivatives.adjustment_w_lin_scratch;
    auto& linear_delta_phi_state = cf.dm_raw;
    if (adjustment_active) {
        MultiFab::Copy(linear_delta_eta, cf.dp, 0, 0, 1, 0);
        applySpectralRieszMultiplier(
            linear_delta_eta,
            ctx.child_operator.fft,
            ctx.child_operator.geomdata,
            1.0,
            1,
            false);
        finalizeCompactDNOField(ctx, linear_delta_eta);
        applyCorrectionTaper(ctx, linear_delta_eta);
        dealiasOnChild(ctx, linear_delta_eta);

        MultiFab::Copy(linear_delta_phi_state, cf.de, 0, 0, 1, 0);
        applyCorrectionTaper(ctx, linear_delta_phi_state);
        dealiasOnChild(ctx, linear_delta_phi_state);
    }

    return {adjustment_active, adjustment_weight};
}

// Step 12a: the parent tendency of stage irk interpolated to the child, I(F0).
void interpolateParentStageTendency(const CorrectionStageContext& ctx, int irk)
{
    auto& level_context = ctx.level_context;
    auto& solution = level_context.solutionState();
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;
    const auto& nested_config = level_context.nestedConfig();

    crescent::nested::rk4::interpolateCoreStageRHSToL1(
        solution.rhs_eta[0], irk,
        level_context.coreLevelGeometry(), level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        nested_config.state_interpolation_order,
        coupling_derivatives.core_interp_patch_scratch,
        coupling_derivatives.rhs_eta_core_on_l1);
    crescent::nested::rk4::interpolateCoreStageRHSToL1(
        solution.rhs_phi[0], irk,
        level_context.coreLevelGeometry(), level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        nested_config.state_interpolation_order,
        coupling_derivatives.core_interp_patch_scratch,
        coupling_derivatives.rhs_phi_core_on_l1);
}

// The pieces of the stage tendency at one point that both RHS assemblies
// use: the correction increment (ta, tb) and the child-assembled background
// tendency F~0 (pressure-free; the applied pressure enters separately).
// While the Dommermuth ramp is active, each is blended from its linear part
// L as L + W (F - L).
struct StagePointTendency
{
    Real delta_eta;
    Real delta_phi;
    Real background_eta;
    Real background_phi;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
StagePointTendency stagePointTendency(
    Real ta, Real tb, Real linear_eta, Real linear_phi_state,
    Real e0, Real e0x, Real e0y, Real p0x, Real p0y,
    Real w0, Real wlin0, Real m0d, Real w0sq,
    int use_adjustment, Real ramp, Real g_grav) noexcept
{
    Real delta_eta = ta;
    Real delta_phi = tb;
    if (use_adjustment) {
        const Real delta_eta_linear = linear_eta;
        const Real delta_phi_linear =
            -g_grav * linear_phi_state;
        delta_eta = delta_eta_linear
            + ramp * (delta_eta - delta_eta_linear);
        delta_phi = delta_phi_linear
            + ramp * (delta_phi - delta_phi_linear);
    }
    const Real f0e_c =
        -(p0x * e0x + p0y * e0y)
        + (1.0 + m0d) * w0;
    const Real f0p_c =
        -g_grav * e0
        - 0.5 * (p0x * p0x
                 + p0y * p0y)
        + 0.5 * (1.0 + m0d) * w0sq;
    const Real background_eta = use_adjustment
        ? wlin0 + ramp * (f0e_c - wlin0)
        : f0e_c;
    const Real background_phi = use_adjustment
        ? -g_grav * e0
            + ramp * (f0p_c + g_grav * e0)
        : f0p_c;
    return {delta_eta, delta_phi, background_eta, background_phi};
}

// Step 12b, the stage tendency. The compact correction tendency
// c = chi * (F~0 - I(F0)) + DeltaF is formed in term_a/term_b, where F~0 is the
// HOS tendency assembled on the child from the interpolated parent fields and
// the parent-dealiased products (pointwise, never transformed); it differs
// from I(F0) by the resolution commutator, which seeds the correction. The
// source part is tapered so that c is compactly supported. c is projected onto
// the added band and added to the interpolated parent tendency:
// RHS = I(F0) + P_add[c] - p1. The child's shared band therefore evolves
// exactly as the interpolated parent; the correction carries only the modes
// the parent cannot represent. I(F0) already contains the parent's applied
// pressure -I(p0), so with a nonzero pressure the source and the RHS use the
// pressure-free I(F0) + I(p0): the child receives -p1 once, after all
// projections, and the correction -(p1 - I(p0)). Until 2026-10-15 the RHS
// subtracted p1 on top of -I(p0) (pressure applied twice); without pressure
// (every standalone run) the arithmetic is unchanged.
void assembleStageRHS(
    const CorrectionStageContext& ctx,
    int irk,
    const StageAdjustment& adjustment)
{
    const int lev = ctx.lev;
    auto& solution = ctx.level_context.solutionState();
    auto& cf = ctx.cf;
    auto& coupling_derivatives =
        ctx.level_context.nestedRuntime().coupling_derivatives;
    auto& linear_delta_eta = coupling_derivatives.adjustment_w_lin_scratch;
    auto& linear_delta_phi_state = cf.dm_raw;
    const bool adjustment_active = adjustment.active;
    const Real adjustment_weight = adjustment.weight;
    const Real g_grav = crescent::constants::GRAVITY;
    amrex::MultiFab rhs_eta_stage(
        solution.rhs_eta[lev], amrex::make_alias, irk, 1);
    amrex::MultiFab rhs_phi_stage(
        solution.rhs_phi[lev], amrex::make_alias, irk, 1);

    // Applied pressure: zero unless an external surface exchange sets it.
    const bool applied_pressure = solution.pressure[0].norm0() > 0.0 ||
                                  solution.pressure[lev].norm0() > 0.0;
    amrex::MultiFab pressure_core_on_l1;
    if (applied_pressure) {
        pressure_core_on_l1.define(cf.term_b.boxArray(),
                                   cf.term_b.DistributionMap(), 1, 0);
        interpolateToChild(ctx.level_context, solution.pressure[0],
                           pressure_core_on_l1);
    }
    const int use_pressure = applied_pressure ? 1 : 0;

    for (MFIter mfi(cf.term_a, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto ta = cf.term_a.array(mfi);
        auto tb = cf.term_b.array(mfi);
        auto f0e = coupling_derivatives.rhs_eta_core_on_l1.const_array(mfi);
        auto f0p = coupling_derivatives.rhs_phi_core_on_l1.const_array(mfi);
        auto ip0 = applied_pressure ? pressure_core_on_l1.const_array(mfi)
                                    : f0p;
        auto s = ctx.chi_source_sigma.const_array(mfi);
        auto e0 = cf.bg_on_l1.eta.const_array(mfi);
        auto e0x = cf.bg_on_l1.eta_x.const_array(mfi);
        auto e0y = cf.bg_on_l1.eta_y.const_array(mfi);
        auto p0x = cf.bg_on_l1.phi_x.const_array(mfi);
        auto p0y = cf.bg_on_l1.phi_y.const_array(mfi);
        auto w0 = cf.bg_on_l1.phi_z.const_array(mfi);
        auto wlin0 = cf.bg_on_l1.w_lin.const_array(mfi);
        auto m0d = cf.bg_on_l1.m0_d.const_array(mfi);
        auto w0sq = cf.bg_on_l1.w0sq_d.const_array(mfi);
        auto linear_de = linear_delta_eta.const_array(mfi);
        auto linear_dp_state = linear_delta_phi_state.const_array(mfi);
        const int use_adjustment = adjustment_active ? 1 : 0;
        const Real ramp = adjustment_weight;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto t = stagePointTendency(
                ta(i, j, k), tb(i, j, k),
                use_adjustment ? linear_de(i, j, k) : 0.0,
                use_adjustment ? linear_dp_state(i, j, k) : 0.0,
                e0(i, j, k), e0x(i, j, k), e0y(i, j, k),
                p0x(i, j, k), p0y(i, j, k), w0(i, j, k), wlin0(i, j, k),
                m0d(i, j, k), w0sq(i, j, k),
                use_adjustment, ramp, g_grav);
            const Real chi = 1.0 - s(i, j, k);
            // pressure-free I(F0_psi): add back the parent's -I(p0)
            const Real f0p_free = use_pressure ? f0p(i, j, k) + ip0(i, j, k)
                                               : f0p(i, j, k);
            ta(i, j, k) = chi * (t.background_eta - f0e(i, j, k)) + t.delta_eta;
            tb(i, j, k) = chi * (t.background_phi - f0p_free) + t.delta_phi;
        });
    }
    applyAddedBandProjector(ctx, cf.term_a);
    applyAddedBandProjector(ctx, cf.term_b);
    for (MFIter mfi(rhs_eta_stage, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto re = rhs_eta_stage.array(mfi);
        auto rp = rhs_phi_stage.array(mfi);
        auto f0e = coupling_derivatives.rhs_eta_core_on_l1.const_array(mfi);
        auto f0p = coupling_derivatives.rhs_phi_core_on_l1.const_array(mfi);
        auto ta = cf.term_a.const_array(mfi);
        auto tb = cf.term_b.const_array(mfi);
        auto p1 = solution.pressure[lev].const_array(mfi);
        auto ip0 = applied_pressure ? pressure_core_on_l1.const_array(mfi)
                                    : f0p;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            re(i, j, k) = f0e(i, j, k) + ta(i, j, k);
            const Real f0p_free = use_pressure ? f0p(i, j, k) + ip0(i, j, k)
                                               : f0p(i, j, k);
            rp(i, j, k) = f0p_free + tb(i, j, k) - p1(i, j, k);
        });
    }
}

} // namespace

void runL1CorrectionFieldStage(
    crescent::shared::LevelContext& level_context,
    int irk,
    amrex::Real stage_time,
    amrex::Real temporal_weight)
{
    const auto ctx = makeCorrectionStageContext(level_context);
    const Real dt = level_context.runtime().dt[ctx.lev];
    const Real fac = temporal_weight * dt;

    formStageCorrection(ctx, irk, fac);
    // Step 2: the stage background (captured post-filter parent fields).
    interpolateStageBackground(ctx, ctx.cf.core_stage[irk]);
    differentiateCorrection(ctx);

    if (ctx.child_operator.perturbation_order <= 1) {
        assembleLinearStageRHS(ctx, irk);
        return;
    }

    // Step 6 (M >= 2): the cross terms, compact support restored, and the
    // parent-parity final projection bracketed by tapers.
    accumulateCrossTerms(ctx);
    finalizeCompactDNOField(ctx, ctx.cf.dw);

    formMetricIncrement(ctx);
    formWSquareIncrement(ctx);
    formMetricWSquareTerm(ctx);
    assembleTendencyIncrement(ctx);
    taperAndDealiasIncrement(ctx);

    const auto adjustment = prepareStageAdjustment(ctx, stage_time);
    interpolateParentStageTendency(ctx, irk);
    assembleStageRHS(ctx, irk, adjustment);
}

void applyPostRK4CorrectionFieldTaper(crescent::shared::LevelContext& level_context)
{
    auto ctx = makeCorrectionStageContext(level_context);
    const int lev = ctx.lev;
    auto& solution = level_context.solutionState();
    auto& cf = ctx.cf;

    // Baseline u0_end = I(parent final state).
    interpolateToChild(level_context, solution.eta[0], cf.cross_a);
    interpolateToChild(level_context, solution.phi[0], cf.cross_b);

    formTaperedCorrection(
        cf.de, solution.eta[lev], cf.cross_a, ctx.chi_source_sigma);
    formTaperedCorrection(
        cf.dp, solution.phi[lev], cf.cross_b, ctx.chi_source_sigma);

    // Remove the shared-band content that the taper and the Runge-Kutta
    // combination can leave in the correction.
    applyAddedBandProjector(ctx, cf.de);
    applyAddedBandProjector(ctx, cf.dp);

    // HL on the compact correction (mean-preserving; no C3, no zero-mean).
    crescent::nested::filter::applyHouLiFilter_L1(
        cf.de,
        cf.dp,
        ctx.child_operator.fft,
        ctx.child_operator.geomdata,
        ctx.child_operator.spectral_settings,
        ctx.child_operator.spectral_workspace.filter_weights[ctx.lev]);

    // u = u0_end + HL[chi * (u - u0_end)]  (single chi; spec §1).
    for (MFIter mfi(solution.eta[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto eta = solution.eta[lev].array(mfi);
        auto phi = solution.phi[lev].array(mfi);
        auto b_eta = cf.cross_a.const_array(mfi);
        auto b_phi = cf.cross_b.const_array(mfi);
        auto de = cf.de.const_array(mfi);
        auto dp = cf.dp.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            eta(i, j, k) = b_eta(i, j, k) + de(i, j, k);
            phi(i, j, k) = b_phi(i, j, k) + dp(i, j, k);
        });
    }
}

void recomputeCorrectionFieldOutputDerivatives(
    crescent::shared::LevelContext& level_context)
{
    auto ctx = makeCorrectionStageContext(level_context);
    const int lev = ctx.lev;
    auto& solution = level_context.solutionState();
    auto& cf = ctx.cf;
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;
    auto& out = coupling_derivatives.l1_blended_derivatives;
    auto& core_on_l1 = coupling_derivatives.core_derivatives_on_l1;

    // Correction after the post-step: u - I(parent final) is the filtered,
    // Fourier-series correction (FFT-safe). Re-taper for exact periodicity.
    interpolateToChild(level_context, solution.eta[0], cf.cross_a);
    interpolateToChild(level_context, solution.phi[0], cf.cross_b);
    formTaperedCorrection(
        cf.de, solution.eta[lev], cf.cross_a, ctx.chi_source_sigma);
    formTaperedCorrection(
        cf.dp, solution.phi[lev], cf.cross_b, ctx.chi_source_sigma);

    auto& workspace = ctx.child_operator.spectral_workspace;
    crescent::physics::computeSurfaceDerivatives(
        solution.eta_x[lev],
        solution.eta_y[lev],
        solution.phi_x[lev],
        solution.phi_y[lev],
        cf.de,
        cf.dp,
        ctx.child_operator.fft,
        ctx.child_operator.geomdata,
        ctx.child_operator.spectral_settings.wavenumbers.dk_x,
        ctx.child_operator.spectral_settings.wavenumbers.dk_y,
        workspace.directional_input_spectral[lev],
        workspace.directional_x_spectral[lev],
        workspace.directional_y_spectral[lev]);
    crescent::dno::computeNestedL1DNO(
        ctx.child_operator,
        cf.de,
        cf.dp,
        solution.wavenumber_pow[lev],
        solution.eta_pow[lev],
        solution.phi_expansion[lev],
        cf.dw,
        false);

    // Order-2 background-correction cross terms with the FINAL parent
    // background, then the same chi re-taper as the stage (spec §4), so the
    // output vertical velocity is the one the stage tendency actually used:
    // w = I(w_L0^{n+1}) + chi * (G_M[de] dp + A1 + A2 + B1 + B2). Until
    // 2026-08-16 the output held only the all-correction DNO G_M[de] dp, so
    // the diagnostic phi_z and the L1 kinetic energy computed from it
    // omitted the cross terms wherever the correction was nonzero.
    const auto& np1 = coupling_derivatives.core_snapshots.np1;
    if (ctx.child_operator.perturbation_order >= 2) {
        const auto core_operator = level_context.coreWaveOperator();
        fillCoreBackgroundSnapshot(
            core_operator,
            solution.eta[0], solution.phi[0],
            np1.eta_x[0], np1.eta_y[0],
            np1.phi_x[0], np1.phi_y[0], np1.phi_z[0],
            cf.core_final);
        interpolateStageBackground(ctx, cf.core_final);
        accumulateCrossTerms(ctx);
    }
    finalizeCompactDNOField(ctx, cf.dw);

    // Output bundle = I(parent np1 bundle) + delta derivatives. The
    // interpolated parent bundle is kept too, so the eta_x_interp and
    // phi_x_interp outputs have their stated meaning.
    interpolateToChild(level_context, np1.eta_x[0], core_on_l1.eta_x[lev]);
    interpolateToChild(level_context, np1.eta_y[0], core_on_l1.eta_y[lev]);
    interpolateToChild(level_context, np1.phi_x[0], core_on_l1.phi_x[lev]);
    interpolateToChild(level_context, np1.phi_y[0], core_on_l1.phi_y[lev]);
    interpolateToChild(level_context, np1.phi_z[0], core_on_l1.phi_z[lev]);
    MultiFab::Copy(out.eta_x[lev], core_on_l1.eta_x[lev], 0, 0, 1, 0);
    MultiFab::Copy(out.eta_y[lev], core_on_l1.eta_y[lev], 0, 0, 1, 0);
    MultiFab::Copy(out.phi_x[lev], core_on_l1.phi_x[lev], 0, 0, 1, 0);
    MultiFab::Copy(out.phi_y[lev], core_on_l1.phi_y[lev], 0, 0, 1, 0);
    MultiFab::Copy(out.phi_z[lev], core_on_l1.phi_z[lev], 0, 0, 1, 0);
    MultiFab::Add(out.eta_x[lev], solution.eta_x[lev], 0, 0, 1, 0);
    MultiFab::Add(out.eta_y[lev], solution.eta_y[lev], 0, 0, 1, 0);
    MultiFab::Add(out.phi_x[lev], solution.phi_x[lev], 0, 0, 1, 0);
    MultiFab::Add(out.phi_y[lev], solution.phi_y[lev], 0, 0, 1, 0);
    MultiFab::Add(out.phi_z[lev], cf.dw, 0, 0, 1, 0);
}

} // namespace crescent::nested::correction_field

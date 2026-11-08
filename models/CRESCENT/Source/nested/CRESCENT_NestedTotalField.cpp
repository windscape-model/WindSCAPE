#include "CRESCENT_NestedTotalField.H"

#include <AMReX_MultiFab.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_Dynamics.H"
#include "CRESCENT_Interpolation.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_NestedC3Dispatch.H"
#include "CRESCENT_NestedDerivativeCoupling.H"
#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_NestedRK4.H"
#include "CRESCENT_NestedSurfaceDerivatives.H"
#include "CRESCENT_NestedTypes.H"
#include "CRESCENT_SolverContext.H"

namespace crescent::nested::total_field {

namespace {

using amrex::MultiFab;

// Copy the live child state, make it periodic with the C3 interface repair,
// and apply the combined child filter. Every child transform consumes such
// a copy; the live state itself is never transformed.
void preparePeriodicCopy(
    crescent::shared::LevelContext& level_context,
    PeriodicDerivativeScratch& periodic_fields)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto mesh = level_context.meshHierarchy();
    auto& nested_runtime = level_context.nestedRuntime();

    MultiFab::Copy(periodic_fields.eta_periodic, solution.eta[lev], 0, 0, 1, 0);
    MultiFab::Copy(periodic_fields.phi_periodic, solution.phi[lev], 0, 0, 1, 0);

    crescent::nested::c3::applyRepairPair(
        lev,
        periodic_fields.eta_periodic,
        periodic_fields.phi_periodic,
        mesh.geom,
        mesh.grids,
        level_context.nestedConfig(),
        nested_runtime.coupling_zones,
        nested_runtime.c3_repair_workspace[lev]);

    const auto level_wave_operator = level_context.waveOperator();
    crescent::nested::filter::applyHouLiFilter_L1(
        periodic_fields.eta_periodic,
        periodic_fields.phi_periodic,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_settings,
        level_wave_operator.spectral_workspace.filter_weights[lev]);
}

void preparePeriodicDerivatives(
    crescent::shared::LevelContext& level_context,
    PeriodicDerivativeScratch& periodic_fields)
{
    preparePeriodicCopy(level_context, periodic_fields);

    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto level_wave_operator = level_context.waveOperator();
    auto& spectral_workspace = level_wave_operator.spectral_workspace;
    auto native_bundle = solution.derivativeBundleAt(lev);

    crescent::nested::derivatives::computeSurfaceDerivatives(
        native_bundle.eta_x,
        native_bundle.eta_y,
        native_bundle.phi_x,
        native_bundle.phi_y,
        periodic_fields.eta_periodic,
        periodic_fields.phi_periodic,
        level_context.levelGeometry(),
        level_wave_operator.fft,
        spectral_workspace.directional_input_spectral[lev],
        spectral_workspace.directional_x_spectral[lev],
        spectral_workspace.directional_y_spectral[lev],
        spectral_workspace.directional_x_valid[lev],
        spectral_workspace.directional_y_valid[lev]);

    crescent::dno::computeNestedL1DNO(
        level_wave_operator,
        periodic_fields.eta_periodic,
        periodic_fields.phi_periodic,
        solution.wavenumber_pow[lev],
        solution.eta_pow[lev],
        solution.phi_expansion[lev],
        native_bundle.phi_z,
        true);
}

} // namespace

void runStage(
    crescent::shared::LevelContext& level_context,
    int irk)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    const auto& nested_config = level_context.nestedConfig();
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;

    auto& periodic_fields =
        level_context.nestedRuntime().periodic_derivative_scratch[lev];

    // Stage 0 evaluates the state of the previous step end. The output
    // preparation of that step (the same copy, C3 repair, filter,
    // derivative and vertical-velocity operations on the same state) is
    // still in the periodic scratch and the native bundle, so it is used
    // instead of a second identical evaluation.
    const bool reuse_step_end_preparation =
        irk == 0 && periodic_fields.holds_step_end_preparation;
    if (!reuse_step_end_preparation) {
        preparePeriodicDerivatives(level_context, periodic_fields);
    }
    periodic_fields.holds_step_end_preparation = false;

    // Tendency blend: form the full nonlinear RHS from the C3-periodized
    // fields (dealiased exactly as on L0), interpolate L0's stage-irk RHS
    // (stage-exact in time), and sigma-blend the two tendencies. All spectral
    // operations act on periodic fields; the blended tendency is non-periodic
    // and is used only in physical space by the RK4 update.
    auto native_l1_bundle =
        level_context.solutionState().derivativeBundleAt(lev);
    const auto& sigma_blend =
        level_context.nestedRuntime().coupling_zones.sigma_blend[lev];
    const auto level_wave_operator = level_context.waveOperator();

    crescent::nested::rk4::computeNativeSurfaceRHS_L1(
        level_wave_operator,
        periodic_fields.eta_periodic,
        native_l1_bundle,
        solution.pressure[lev],
        coupling_derivatives.rhs_eta_native,
        coupling_derivatives.rhs_phi_native);

    if (level_context.core().wave.spectral.adjustment_periods > 0.0) {
        // Same Dommermuth ramp as the core, applied to the native child
        // tendency from the periodized pair (transform-admissible), so
        // that (1-sigma) F_L1 + sigma I(F_L0) stays ramp-consistent. The
        // native RHS fabs are single-component, hence rk_stage 0 here.
        constexpr amrex::Real stage_c[4] = {0.0, 0.5, 0.5, 1.0};
        const auto& wave = level_context.core().wave;
        const amrex::Real t_peak = 2.0 * crescent::constants::PI /
            wave.initialization.peak_frequency;
        crescent::dynamics::applyAdjustmentRampToStageRHS(
            coupling_derivatives.rhs_eta_native,
            coupling_derivatives.rhs_phi_native,
            periodic_fields.eta_periodic,
            periodic_fields.phi_periodic,
            solution.pressure[lev],
            level_wave_operator.fft,
            level_wave_operator.geomdata,
            runtime.t_new[lev] + stage_c[irk] * runtime.dt[lev],
            wave.spectral.adjustment_periods * t_peak,
            0,
            coupling_derivatives.adjustment_w_lin_scratch);
    }

    crescent::nested::rk4::interpolateCoreStageRHSToL1(
        solution.rhs_eta[0],
        irk,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        nested_config.derivative_interpolation_order,
        coupling_derivatives.core_interp_patch_scratch,
        coupling_derivatives.rhs_eta_core_on_l1);
    crescent::nested::rk4::interpolateCoreStageRHSToL1(
        solution.rhs_phi[0],
        irk,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        nested_config.derivative_interpolation_order,
        coupling_derivatives.core_interp_patch_scratch,
        coupling_derivatives.rhs_phi_core_on_l1);

    amrex::MultiFab rhs_eta_stage(
        solution.rhs_eta[lev], amrex::make_alias, irk, 1);
    amrex::MultiFab rhs_phi_stage(
        solution.rhs_phi[lev], amrex::make_alias, irk, 1);
    crescent::nested::rk4::blendWithParent(
        coupling_derivatives.rhs_eta_native,
        coupling_derivatives.rhs_eta_core_on_l1,
        sigma_blend,
        rhs_eta_stage);
    crescent::nested::rk4::blendWithParent(
        coupling_derivatives.rhs_phi_native,
        coupling_derivatives.rhs_phi_core_on_l1,
        sigma_blend,
        rhs_phi_stage);
}


void prepareOutputDerivatives(
    crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& periodic_fields =
        level_context.nestedRuntime().periodic_derivative_scratch[lev];
    preparePeriodicDerivatives(level_context, periodic_fields);
    periodic_fields.holds_step_end_preparation = true;
    crescent::nested::derivative_coupling::prepareNestedL1FinalOutputDerivatives(
        level_context,
        level_context.solutionState().derivativeBundleAt(lev));
}

void applyStateRelaxation(
    crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    auto& nested_runtime = level_context.nestedRuntime();
    // The stage scratch is free here: the next reader is the end-of-step
    // output preparation, which rebuilds it from the relaxed state.
    auto& periodic_fields = nested_runtime.periodic_derivative_scratch[lev];
    auto& post_rk4_scratch = nested_runtime.post_rk4_state_scratch[lev];

    preparePeriodicCopy(level_context, periodic_fields);

    crescent::interpolation::interpolateParentToChild(
        solution.eta[0],
        solution.phi[0],
        post_rk4_scratch.interpolated_core_eta,
        post_rk4_scratch.interpolated_core_phi,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        post_rk4_scratch.coarseInterpPatch(),
        level_context.nestedConfig().state_interpolation_order);

    const MultiFab& sigma_blend = nested_runtime.coupling_zones.sigma_blend[lev];
    for (amrex::MFIter mfi(solution.eta[lev]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto eta_arr = solution.eta[lev].array(mfi);
        auto phi_arr = solution.phi[lev].array(mfi);
        auto eta_filt_arr = periodic_fields.eta_periodic.const_array(mfi);
        auto phi_filt_arr = periodic_fields.phi_periodic.const_array(mfi);
        auto eta_core_arr = post_rk4_scratch.interpolated_core_eta.const_array(mfi);
        auto phi_core_arr = post_rk4_scratch.interpolated_core_phi.const_array(mfi);
        auto sigma_arr = sigma_blend.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real sig = sigma_arr(i, j, k);
            eta_arr(i, j, k) =
                (1.0 - sig) * eta_filt_arr(i, j, k) + sig * eta_core_arr(i, j, k);
            phi_arr(i, j, k) =
                (1.0 - sig) * phi_filt_arr(i, j, k) + sig * phi_core_arr(i, j, k);
        });
    }
}

} // namespace crescent::nested::total_field

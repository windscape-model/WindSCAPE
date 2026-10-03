#include "CRESCENT_NestedRK4Stage.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_Dynamics.H"
#include "CRESCENT_NestedDerivativeBundles.H"
#include "CRESCENT_NestedDerivativeCoupling.H"
#include "CRESCENT_NestedPeriodicSurfaceDerivatives.H"
#include "CRESCENT_NestedRuntimeLogging.H"
#include "CRESCENT_NestedStageDiagnostics.H"
#include "CRESCENT_NestedStageResearchDiagnostics.H"
#include "CRESCENT_NestedTypes.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"

namespace {

using namespace crescent::constants;
namespace derivative_coupling = crescent::nested::derivative_coupling;
namespace nested_stage_diagnostics = crescent::nested::stage_diagnostics;
namespace nested_stage_research = crescent::nested::research_diagnostics;
namespace runtime_logging = crescent::nested::runtime_logging;

} // namespace

namespace crescent::nested::rk4_stage {

void runL1RK4Stage(
    crescent::shared::LevelEvolutionContext& level_context,
    int irk,
    amrex::Real stage_time,
    amrex::Real temporal_weight)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    const auto& nested_config = level_context.nestedConfig();
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;

    crescent::dynamics::advanceSurfaceRK4(
        solution.eta[lev],
        solution.phi[lev],
        solution.eta_saved[lev],
        solution.phi_saved[lev],
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        runtime.dt[lev],
        irk);

    auto stage_diagnostics_state =
        nested_stage_diagnostics::makeRK4StageDiagnosticsContext(
            level_context,
            irk,
            stage_time,
            temporal_weight);

    const bool should_log_stage_research =
        nested_stage_research::shouldLogDetailedStageDiagnostics(
            stage_diagnostics_state.verbosity);
    const auto& accurate_region =
        level_context.nestedRuntime().coupling_zones.domains[lev].box_A;

    runtime_logging::logNestedL1StageHeader(
        stage_diagnostics_state.irk,
        stage_diagnostics_state.temporal_weight,
        stage_diagnostics_state.verbosity);

    if (should_log_stage_research) {
        nested_stage_research::logNestedL1StageEntryDiagnostics(
            solution.eta[lev],
            accurate_region,
            stage_diagnostics_state.irk,
            stage_diagnostics_state.verbosity);
    }

    auto& periodic_fields =
        level_context.nestedRuntime().periodic_derivative_scratch[lev];

    crescent::nested::periodic_surface_derivatives::prepareRK4StagePeriodicL1Derivatives(
        level_context,
        periodic_fields,
        stage_diagnostics_state);

    nested_stage_diagnostics::reportPreparedRK4Stage(
        stage_diagnostics_state,
        level_context,
        solution.phi[lev],
        periodic_fields.phi_periodic);

    if (should_log_stage_research && stage_diagnostics_state.irk == 0) {
        nested_stage_research::logNestedL1CorePhiZReferenceDiagnostics(
            level_context.coreLevelGeometry(),
            coupling_derivatives.core_snapshots.n.phi_z[0],
            level_context.levelGeometry(),
            accurate_region);
    }

    const auto stage_views =
        crescent::nested::bundles::makeL1StageDerivativeViews(
            coupling_derivatives,
            lev);

    auto native_l1_bundle =
        level_context.solutionState().derivativeBundleAt(lev);
    derivative_coupling::interpolateAndBlendCoreDerivativesToL1(
        stage_views.core_n,
        stage_views.core_np1,
        temporal_weight,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        nested_config.derivative_interpolation_order,
        native_l1_bundle,
        level_context.nestedRuntime().coupling_zones.sigma_blend[lev],
        nested_config.blending_correction_y_smoothing_passes,
        coupling_derivatives.core_stage_derivative_scratch,
        coupling_derivatives.core_interp_patch_scratch,
        coupling_derivatives.blend_correction_delta_raw,
        coupling_derivatives.blend_correction_delta_smooth,
        stage_views.core_derivatives_on_l1,
        stage_views.blended_derivatives);

    if (should_log_stage_research) {
        nested_stage_research::logNestedL1InterpolatedDerivativeDiagnostics(
            stage_views.core_derivatives_on_l1,
            accurate_region);
        nested_stage_research::logNestedL1BlendedDerivativeDiagnostics(
            stage_views.blended_derivatives);
    }

    crescent::nested::rk4::computeSurfaceRHS_L1(
        solution.rhs_eta[lev],
        solution.rhs_phi[lev],
        solution.eta[lev],
        stage_views.blended_derivatives.eta_x,
        stage_views.blended_derivatives.eta_y,
        stage_views.blended_derivatives.phi_x,
        stage_views.blended_derivatives.phi_y,
        stage_views.blended_derivatives.phi_z,
        irk,
        GRAVITY);

    if (should_log_stage_research) {
        nested_stage_research::logNestedL1RhsDiagnostics(
            solution.rhs_eta[lev],
            solution.rhs_phi[lev],
            accurate_region,
            stage_diagnostics_state.irk);
    }
}

} // namespace crescent::nested::rk4_stage

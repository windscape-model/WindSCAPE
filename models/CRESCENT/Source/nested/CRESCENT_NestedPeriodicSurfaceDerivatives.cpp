// Shared periodic-scratch preparation for nested RK4 stages and output fields.

#include "CRESCENT_NestedPeriodicSurfaceDerivatives.H"

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_NestedC3Dispatch.H"
#include "CRESCENT_NestedDerivativeBundles.H"
#include "CRESCENT_NestedDerivativeExecution.H"
#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_NestedPeriodicResearchDiagnostics.H"
#include "CRESCENT_NestedStageDiagnostics.H"
#include "CRESCENT_NestedTypes.H"
#include "CRESCENT_Spectral.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"

namespace {

using LevelMeshView = crescent::shared::LevelMeshView;
using PeriodicDerivativeScratch = crescent::nested::PeriodicDerivativeScratch;
using RK4StageDiagnosticsContext =
    crescent::nested::stage_diagnostics::RK4StageDiagnosticsContext;
using LevelEvolutionContext = crescent::shared::LevelEvolutionContext;
using DerivativeBundleView = crescent::shared::DerivativeBundleView;
using FreeSurfaceStateView = crescent::shared::FreeSurfaceStateView;
using ConstFreeSurfaceStateView = crescent::shared::ConstFreeSurfaceStateView;
using LevelWaveOperatorView = crescent::shared::LevelWaveOperatorView;
using NestedAlgorithmConfig = crescent::nested::AlgorithmConfig;
using CouplingZoneState = crescent::nested::CouplingZoneState;
using amrex::MultiFab;

struct NestedPeriodicDerivativeContext
{
    int lev;
    int istep;
    int verbosity;
    int use_c3_repair;
    int c3_use_interface;
    int c3_repair_after_filter;
    int enable_initial_c3_effect_diagnostic;
    int filter_periodic_copies_before_linear_derivatives;
    int linear_derivative_spectral_filter_mode;

    crescent::shared::MeshHierarchyView mesh;
    LevelMeshView level_mesh;
    LevelWaveOperatorView level_wave_operator;
    const NestedAlgorithmConfig& nested_config;
    const CouplingZoneState& coupling_zones;
    crescent::nested::C3RepairWorkspace& c3_repair_workspace;

    ConstFreeSurfaceStateView source_state;
    const amrex::Box& spectral_extract_box;

    DerivativeBundleView derivative_bundle;
    const MultiFab& wavenumber_pow;
    MultiFab& eta_pow;
    MultiFab& phi_expansion;
};

struct PeriodicScratchPreparationPlan
{
    bool allow_initial_c3_effect_diagnostic;
    bool require_interface_c3_repair;
    int log_stage;
    int dno_verbosity;
    bool cleanup_stage_phi_z;
    MultiFab* phi_post_c3_snapshot;
};

NestedPeriodicDerivativeContext makeNestedPeriodicDerivativeContext(
    LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());
    const auto& nested_config = level_context.nestedConfig();
    auto& coupling_zones = level_context.nestedRuntime().coupling_zones;

    return {
        lev,
        level_context.runtime().istep[lev],
        level_context.runtime().verbosity,
        nested_config.use_c3_repair,
        nested_config.c3_use_interface,
        nested_config.c3_repair_after_filter,
        nested_config.enable_initial_c3_effect_diagnostic,
        nested_config.filter_periodic_copies_before_linear_derivatives,
        nested_config.linear_derivative_spectral_filter_mode,
        level_context.meshHierarchy(),
        level_context.levelMesh(),
        level_wave_operator,
        nested_config,
        coupling_zones,
        level_context.nestedRuntime().c3_repair_workspace[lev],
        crescent::shared::asConstFreeSurfaceState(
            level_context.solutionState().freeSurfaceStateAt(lev)),
        coupling_zones.domains[lev].box_S,
        level_context.solutionState().derivativeBundleAt(lev),
        solution.wavenumber_pow[lev],
        solution.eta_pow[lev],
        solution.phi_expansion[lev]
    };
}

PeriodicScratchPreparationPlan makeRK4StagePeriodicScratchPreparationPlan(
    RK4StageDiagnosticsContext& stage_diagnostics)
{
    return {
        true,
        false,
        stage_diagnostics.irk,
        stage_diagnostics.verbosity,
        true,
        stage_diagnostics.phiStageSnapshot()
    };
}

PeriodicScratchPreparationPlan makeFinalOutputPeriodicScratchPreparationPlan()
{
    return {
        false,
        true,
        -1,
        0,
        false,
        nullptr
    };
}

void applyNestedPeriodicC3Repair(
    const NestedPeriodicDerivativeContext& context,
    FreeSurfaceStateView state,
    int repair_verbosity)
{
    crescent::nested::c3::applyRepairPair(
        context.lev,
        state.eta,
        state.phi,
        context.mesh.geom,
        context.mesh.grids,
        context.nested_config,
        context.coupling_zones,
        repair_verbosity,
        context.c3_repair_workspace);
}

void copyFreeSurfaceState(
    ConstFreeSurfaceStateView source_state,
    FreeSurfaceStateView target_state)
{
    amrex::MultiFab::Copy(target_state.eta, source_state.eta, 0, 0, 1, 0);
    amrex::MultiFab::Copy(target_state.phi, source_state.phi, 0, 0, 1, 0);
}

void capturePeriodicPostC3SnapshotIfRequested(
    const PeriodicScratchPreparationPlan& plan,
    const PeriodicDerivativeScratch& periodic_fields)
{
    if (plan.phi_post_c3_snapshot == nullptr) {
        return;
    }

    // Capture the post-C3 state before periodic filtering mutates the stage copy.
    amrex::MultiFab::Copy(
        *plan.phi_post_c3_snapshot,
        periodic_fields.phi_periodic,
        0,
        0,
        1,
        0);
}

void copyCurrentNestedL1StateToPeriodicScratch(
    const NestedPeriodicDerivativeContext& context,
    PeriodicDerivativeScratch& periodic_fields,
    const PeriodicScratchPreparationPlan& plan)
{
    const bool run_initial_c3_effect_diagnostic =
        plan.allow_initial_c3_effect_diagnostic &&
        context.use_c3_repair == 1 &&
        context.istep == 0 &&
        context.verbosity >= 2 &&
        context.enable_initial_c3_effect_diagnostic == 1;

    if (run_initial_c3_effect_diagnostic) {
        crescent::nested::research_diagnostics::runInitialC3RepairEffectDiagnostic(
            context.level_mesh,
            context.source_state.eta,
            context.source_state.phi,
            periodic_fields,
            [&](amrex::MultiFab& eta_mf, amrex::MultiFab& phi_mf, int repair_verbosity) {
                applyNestedPeriodicC3Repair(
                    context,
                    {eta_mf, phi_mf},
                    repair_verbosity);
            });
        return;
    }

    copyFreeSurfaceState(
        context.source_state,
        periodic_fields.periodicState());

    const bool should_apply_periodic_repair =
        context.use_c3_repair == 1 &&
        !(plan.require_interface_c3_repair && context.c3_use_interface != 1);
    if (!should_apply_periodic_repair) {
        return;
    }

    if (context.verbosity > 1) {
        crescent::nested::research_diagnostics::logNestedPreRepairYVariation(
            context.source_state.eta,
            context.verbosity);
    }

    applyNestedPeriodicC3Repair(
        context,
        periodic_fields.periodicState(),
        0);

    if (context.verbosity > 1) {
        crescent::nested::research_diagnostics::logNestedPostRepairYVariation(
            periodic_fields.eta_periodic,
            context.verbosity);
    }
}

void filterNestedL1PeriodicScratch(
    const NestedPeriodicDerivativeContext& context,
    PeriodicDerivativeScratch& periodic_fields,
    const PeriodicScratchPreparationPlan& plan)
{
    periodic_fields.snapshotLinearPrefilterCopies();

    if (context.verbosity >= 2 && plan.log_stage == 0) {
        crescent::nested::research_diagnostics::logNestedPeriodicPrefilterState(
            periodic_fields.eta_periodic,
            plan.log_stage,
            context.verbosity);
    }

    crescent::nested::filter::applyHouLiFilter_L1(
        periodic_fields.eta_periodic,
        periodic_fields.phi_periodic,
        context.level_wave_operator.fft,
        context.level_wave_operator.geomdata,
        context.level_wave_operator.spectral_settings,
        context.verbosity >= 2 && context.istep == 0 && plan.log_stage == 0);

    if (context.use_c3_repair == 1 && context.c3_repair_after_filter == 1) {
        applyNestedPeriodicC3Repair(
            context,
            periodic_fields.periodicState(),
            0);
    }

    if (context.verbosity >= 2 && plan.log_stage == 0) {
        crescent::nested::research_diagnostics::logNestedPeriodicPostfilterState(
            periodic_fields.eta_periodic,
            plan.log_stage,
            context.verbosity);
    }
}

void cleanupNestedL1StagePhiZ(
    const NestedPeriodicDerivativeContext& context,
    MultiFab& phi_z_filter_scratch)
{
    crescent::nested::filter::applyHouLiFilter_L1(
        context.derivative_bundle.phi_z,
        phi_z_filter_scratch,
        context.level_wave_operator.fft,
        context.level_wave_operator.geomdata,
        context.level_wave_operator.spectral_settings,
        false);
    crescent::spectral::dealiasing(
        context.derivative_bundle.phi_z,
        context.level_wave_operator.fft,
        context.level_wave_operator.geomdata,
        context.level_wave_operator.spectral_workspace.dealiasing_scratch[context.level_wave_operator.lev]);
}

void computePreparedNestedL1PeriodicDerivatives(
    const NestedPeriodicDerivativeContext& context,
    PeriodicDerivativeScratch& periodic_fields,
    const PeriodicScratchPreparationPlan& plan)
{
    const auto linear_derivative_inputs =
        crescent::nested::bundles::selectLinearDerivativeSources(
            periodic_fields.view(),
            context.filter_periodic_copies_before_linear_derivatives == 1);
    const auto periodic_source_state =
        crescent::shared::asConstFreeSurfaceState(periodic_fields.periodicState());

    crescent::nested::derivative_execution::computePreparedPeriodicDerivativeBundle(
        {
            context.level_wave_operator,
            context.level_mesh.geom,
            linear_derivative_inputs,
            periodic_source_state,
            context.spectral_extract_box,
            context.linear_derivative_spectral_filter_mode,
            context.verbosity,
            plan.dno_verbosity,
            context.derivative_bundle,
            context.wavenumber_pow,
            context.eta_pow,
            context.phi_expansion
        });

    if (!plan.cleanup_stage_phi_z) {
        return;
    }

    cleanupNestedL1StagePhiZ(context, periodic_fields.phi_z_filter_scratch);

    if (context.verbosity >= 2) {
        crescent::nested::research_diagnostics::logNestedL1PhiZCleanupDiagnostics(
            context.derivative_bundle.eta_x,
            context.derivative_bundle.phi_z);
    }
}

void preparePeriodicScratchDerivatives(
    const NestedPeriodicDerivativeContext& context,
    PeriodicDerivativeScratch& periodic_fields,
    const PeriodicScratchPreparationPlan& plan)
{
    copyCurrentNestedL1StateToPeriodicScratch(
        context,
        periodic_fields,
        plan);

    capturePeriodicPostC3SnapshotIfRequested(plan, periodic_fields);

    filterNestedL1PeriodicScratch(
        context,
        periodic_fields,
        plan);

    computePreparedNestedL1PeriodicDerivatives(
        context,
        periodic_fields,
        plan);
}

} // namespace

namespace crescent::nested::periodic_surface_derivatives {

void prepareRK4StagePeriodicL1Derivatives(
    LevelEvolutionContext& level_context,
    PeriodicDerivativeScratch& periodic_fields,
    RK4StageDiagnosticsContext& stage_diagnostics)
{
    const auto context =
        makeNestedPeriodicDerivativeContext(level_context);
    const auto preparation_plan =
        makeRK4StagePeriodicScratchPreparationPlan(stage_diagnostics);

    preparePeriodicScratchDerivatives(
        context,
        periodic_fields,
        preparation_plan);
}

void prepareFinalOutputPeriodicL1Derivatives(
    LevelEvolutionContext& level_context,
    PeriodicDerivativeScratch& periodic_fields)
{
    const auto context =
        makeNestedPeriodicDerivativeContext(level_context);
    const auto preparation_plan =
        makeFinalOutputPeriodicScratchPreparationPlan();

    preparePeriodicScratchDerivatives(
        context,
        periodic_fields,
        preparation_plan);
}

} // namespace crescent::nested::periodic_surface_derivatives

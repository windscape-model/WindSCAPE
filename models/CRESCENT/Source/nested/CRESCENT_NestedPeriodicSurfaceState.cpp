// Shared post-RK4 periodic-state preparation and core-state blending.

#include "CRESCENT_NestedPeriodicSurfaceState.H"

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_Interpolation.H"
#include "CRESCENT_NestedC3Dispatch.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_NestedRuntimeLogging.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"
#include "shared/CRESCENT_LevelWaveOperatorView.H"

namespace {

using FreeSurfaceStateView = crescent::shared::FreeSurfaceStateView;
using ConstFreeSurfaceStateView = crescent::shared::ConstFreeSurfaceStateView;
using LevelWaveSpectralSettings = crescent::shared::LevelWaveSpectralSettings;
using LevelMeshView = crescent::shared::LevelMeshView;
using MeshHierarchyView = crescent::shared::MeshHierarchyView;
using LevelEvolutionContext = crescent::shared::LevelEvolutionContext;
using CouplingZoneState = crescent::nested::CouplingZoneState;
using NestedAlgorithmConfig = crescent::nested::AlgorithmConfig;

struct NestedPeriodicStatePreparationContext
{
    int lev;
    int verbosity;
    int istep;
    int c3_repair_after_filter;
    LevelWaveSpectralSettings spectral_settings;
    MeshHierarchyView mesh;
    LevelMeshView level_mesh;
    const NestedAlgorithmConfig& nested_config;
    const CouplingZoneState& coupling_zones;
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft;
    ConstFreeSurfaceStateView source_state;
    crescent::nested::C3RepairWorkspace& c3_repair_workspace;
};

struct NestedCoreStateInterpolationContext
{
    int interpolation_order;
    const amrex::Geometry& core_geom;
    const amrex::Geometry& nested_geom;
    amrex::IntVect core_to_l1_ref_ratio;
    ConstFreeSurfaceStateView core_state;
};

struct NestedStateBlendTargetContext
{
    const amrex::MultiFab& sigma_blend;
    FreeSurfaceStateView output_state;
};

NestedPeriodicStatePreparationContext makeNestedPeriodicStatePreparationContext(
    LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());

    return {
        lev,
        level_context.runtime().verbosity,
        level_context.runtime().istep[lev],
        level_context.nestedConfig().c3_repair_after_filter,
        level_wave_operator.spectral_settings,
        level_context.meshHierarchy(),
        level_context.levelMesh(),
        level_context.nestedConfig(),
        level_context.nestedRuntime().coupling_zones,
        level_wave_operator.fft,
        crescent::shared::asConstFreeSurfaceState(
            level_context.solutionState().freeSurfaceStateAt(lev)),
        level_context.nestedRuntime().c3_repair_workspace[lev]
    };
}

NestedCoreStateInterpolationContext makeNestedCoreStateInterpolationContext(
    LevelEvolutionContext& level_context)
{
    return {
        level_context.nestedConfig().state_interpolation_order,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        crescent::shared::asConstFreeSurfaceState(
            level_context.solutionState().freeSurfaceStateAt(0))
    };
}

NestedStateBlendTargetContext makeNestedStateBlendTargetContext(
    LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    return {
        level_context.nestedRuntime().coupling_zones.sigma_blend[lev],
        level_context.solutionState().freeSurfaceStateAt(lev)
    };
}

void copyFreeSurfaceState(
    ConstFreeSurfaceStateView source_state,
    FreeSurfaceStateView target_state)
{
    const int nghost = target_state.eta.nGrow();
    amrex::MultiFab::Copy(target_state.eta, source_state.eta, 0, 0, 1, nghost);
    amrex::MultiFab::Copy(target_state.phi, source_state.phi, 0, 0, 1, nghost);
}

void preparePeriodicL1StateScratch(
    const NestedPeriodicStatePreparationContext& preparation_context,
    FreeSurfaceStateView periodic_state)
{
    copyFreeSurfaceState(
        preparation_context.source_state,
        periodic_state);

    crescent::nested::runtime_logging::logNestedStateRangePair(
        "Before C³-repair",
        periodic_state.eta,
        periodic_state.phi,
        preparation_context.verbosity);

    crescent::nested::c3::applyRepairPair(
        preparation_context.lev,
        periodic_state.eta,
        periodic_state.phi,
        preparation_context.mesh.geom,
        preparation_context.mesh.grids,
        preparation_context.nested_config,
        preparation_context.coupling_zones,
        0,
        preparation_context.c3_repair_workspace);

    crescent::nested::runtime_logging::logNestedStateRangePair(
        "After C³-repair",
        periodic_state.eta,
        periodic_state.phi,
        preparation_context.verbosity);

    crescent::nested::filter::applyHouLiFilter_L1(
        periodic_state.eta,
        periodic_state.phi,
        preparation_context.fft,
        preparation_context.level_mesh.domain(),
        preparation_context.spectral_settings,
        preparation_context.verbosity >= 2 &&
            preparation_context.istep == 0);

    if (preparation_context.nested_config.use_c3_repair == 1 &&
        preparation_context.c3_repair_after_filter == 1) {
        crescent::nested::c3::applyRepairPair(
            preparation_context.lev,
            periodic_state.eta,
            periodic_state.phi,
            preparation_context.mesh.geom,
            preparation_context.mesh.grids,
            preparation_context.nested_config,
            preparation_context.coupling_zones,
            0,
            preparation_context.c3_repair_workspace);
    }

    crescent::nested::runtime_logging::logNestedStateRangePair(
        "After filter",
        periodic_state.eta,
        periodic_state.phi,
        preparation_context.verbosity);
}

void interpolateCoreStateToL1(
    const NestedCoreStateInterpolationContext& interpolation_context,
    FreeSurfaceStateView core_state_on_l1,
    amrex::MultiFab& coarse_interp_patch)
{
    crescent::interpolation::interpolateParentToChild(
        interpolation_context.core_state,
        core_state_on_l1,
        interpolation_context.core_geom,
        interpolation_context.nested_geom,
        interpolation_context.core_to_l1_ref_ratio,
        coarse_interp_patch,
        interpolation_context.interpolation_order);
}

void blendFilteredL1StateWithCore(
    const NestedStateBlendTargetContext& blend_target_context,
    ConstFreeSurfaceStateView filtered_l1_state,
    ConstFreeSurfaceStateView interpolated_core_state)
{
    for (amrex::MFIter mfi(blend_target_context.output_state.eta);
         mfi.isValid();
         ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto eta_arr = blend_target_context.output_state.eta.array(mfi);
        auto phi_arr = blend_target_context.output_state.phi.array(mfi);
        auto eta_filt_arr = filtered_l1_state.eta.const_array(mfi);
        auto phi_filt_arr = filtered_l1_state.phi.const_array(mfi);
        auto eta_core_arr = interpolated_core_state.eta.const_array(mfi);
        auto phi_core_arr = interpolated_core_state.phi.const_array(mfi);
        auto sigma_arr = blend_target_context.sigma_blend.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real sig = sigma_arr(i, j, k);
            eta_arr(i, j, k) =
                (1.0 - sig) * eta_filt_arr(i, j, k) + sig * eta_core_arr(i, j, k);
            phi_arr(i, j, k) =
                (1.0 - sig) * phi_filt_arr(i, j, k) + sig * phi_core_arr(i, j, k);
        });
    }
}

} // namespace

namespace crescent::nested::periodic_surface_state {

void applyPostRK4CoreStateBlend(
    LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    auto& solution = level_context.solutionState();
    const auto& runtime = level_context.runtime();
    const auto periodic_state_preparation_context =
        makeNestedPeriodicStatePreparationContext(level_context);
    const auto core_state_interpolation_context =
        makeNestedCoreStateInterpolationContext(level_context);
    const auto state_blend_target_context =
        makeNestedStateBlendTargetContext(level_context);

    crescent::nested::runtime_logging::logNestedStateFilterBlendStart(
        runtime.verbosity);

    auto& post_rk4_scratch =
        level_context.nestedRuntime().post_rk4_state_scratch[lev];
    auto periodic_l1_state = post_rk4_scratch.periodicL1View();
    auto interpolated_core_state = post_rk4_scratch.interpolatedCoreView();

    preparePeriodicL1StateScratch(
        periodic_state_preparation_context,
        periodic_l1_state);

    interpolateCoreStateToL1(
        core_state_interpolation_context,
        interpolated_core_state,
        post_rk4_scratch.coarseInterpPatch());

    crescent::nested::runtime_logging::logNestedStateRangePair(
        "L0 interpolated",
        interpolated_core_state.eta,
        interpolated_core_state.phi,
        runtime.verbosity);

    blendFilteredL1StateWithCore(
        state_blend_target_context,
        post_rk4_scratch.periodicL1ConstView(),
        post_rk4_scratch.interpolatedCoreConstView());

    crescent::nested::runtime_logging::logNestedStateFilterBlendComplete(
        solution.eta[lev],
        solution.phi[lev],
        runtime.verbosity);
}

} // namespace crescent::nested::periodic_surface_state

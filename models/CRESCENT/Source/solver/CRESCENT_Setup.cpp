#include "CRESCENT_Setup.H"

#include <algorithm>
#include <cmath>
#include <string>

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "CRESCENT_CoreWaveInitialization.H"
#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedCoreSnapshots.H"
#include "CRESCENT_NestedSetup.H"
#include "CRESCENT_NestedTimestep.H"
#include "CRESCENT_Physics.H"
#include "CRESCENT_PostStep.H"
#include "CRESCENT_SetupParameters.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_Spectral.H"
#include "CRESCENT_SpectralGridMetrics.H"

namespace crescent::solver::setup {

namespace {

void computeLevelWavenumberSpacing(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& wave = level_context.wave();
    const int verbosity = level_context.runtime().verbosity;
    const auto level_mesh = level_context.levelMesh();
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(
            level_mesh.grids.minimalBox(),
            level_mesh.geom);
    auto& wavenumbers = wave.levelWavenumbers(lev);
    wavenumbers = spectral_grid.wavenumberSpacing();

    if (verbosity >= 1) {
        amrex::Print() << "Computed wavenumber spacing: dk_x[" << lev << "] = "
                       << wavenumbers.dk_x
                       << ", dk_y[" << lev << "] = "
                       << wavenumbers.dk_y
                       << " (Lx=" << spectral_grid.Lx
                       << ", Ly=" << spectral_grid.Ly << ")" << std::endl;
    }
}

void computeLevelWavenumberPowers(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    const auto& wave = level_context.wave();
    const auto& wavenumbers = wave.levelWavenumbers(lev);
    crescent::spectral::computeWavenumberPowers(
        core.solution.wavenumber_pow[lev],
        wavenumbers.dk_x,
        wavenumbers.dk_y);
}

void initializeLevelSpectralMetrics(
    const crescent::shared::LevelContext& level_context)
{
    computeLevelWavenumberSpacing(level_context);
    computeLevelWavenumberPowers(level_context);
}

void prepareLevelWaveState(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    computeLevelWavenumberSpacing(level_context);

    if (lev == 0) {
        crescent::init::initializeCoreWaveField(
            core,
            level_context.levelMesh(),
            level_context.finestLevel());
    } else {
        AMREX_ALWAYS_ASSERT(lev > 0);
        if (level_context.runtime().verbosity >= 1) {
            amrex::Print() << "\n=== Skipping wave initialization for L" << lev
                           << " (will be interpolated from parent) ===" << std::endl;
        }
    }

    computeLevelWavenumberPowers(level_context);
}

void applyInitializationWaveFilter(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    if (lev > 0) {
        // The child state is still the zero-initialized storage; it is
        // interpolated from the filtered parent afterwards.
        return;
    }
    auto& core = level_context.core();
    const auto level_wave_operator = level_context.waveOperator();
    crescent::spectral::applyLowPassFilter(
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        core.solution.eta[lev],
        core.solution.phi[lev],
        level_wave_operator.spectral_settings,
        level_wave_operator.spectral_workspace.filter_weights[lev]);
}

void initializeCoreDerivedFields(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev == 0);
    auto& core = level_context.core();
    const auto level_wave_operator = level_context.waveOperator();
    auto derivative_bundle = core.solution.derivativeBundleAt(lev);
    auto& solution = core.solution;

    crescent::dno::computeSolutionStateDerivativeBundle(
        level_wave_operator,
        core,
        derivative_bundle);

    crescent::physics::computeHorizontalVelocities(
        solution.u[lev],
        solution.v[lev],
        derivative_bundle,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev]);
}

void initializeLevelSpectralWorkspace(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    auto& root = core.root_layouts;
    auto& workspace = core.spectral_workspace;
    const auto& wave = level_context.wave();
    const auto level_mesh = level_context.levelMesh();

    amrex::Box fft_box = level_mesh.grids.minimalBox();
    workspace.fft[lev] =
        std::make_unique<amrex::FFT::R2C<amrex::Real>>(fft_box);

    auto const& [cba, cdm] = workspace.fft[lev]->getSpectralDataLayout();

    workspace.phi_expansion_spectrum[lev].define(
        cba, cdm, wave.spectral.perturbation_order, 0);

    workspace.compute_rhs_metric[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.compute_rhs_phi_z_sq[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.expansion_accum[lev].define(
        level_mesh.grids, level_mesh.dmap, wave.spectral.perturbation_order, 0);
    workspace.spectral_work[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.expansion_deriv[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.deriv_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_input_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_x_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_y_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_x_valid[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.directional_y_valid[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    // scalar_nghost follows the convention in SolutionState::defineScratchLevelStorage
    // so the scratch covers the widest nGrow that dealiasing consumers pass.
    const int dealiasing_nghost = (lev == 0) ? 1 : 2;
    workspace.dealiasing_scratch[lev].define(
        level_mesh.grids, level_mesh.dmap, 1, dealiasing_nghost);

    core.solution.wavenumber_pow[lev].clear();
    core.solution.wavenumber_pow[lev].define(cba, cdm, wave.spectral.perturbation_order, 0);

    amrex::Vector<int> croot_map(1);
    croot_map[0] = amrex::ParallelDescriptor::IOProcessorNumber();
    root.cba[lev] = amrex::BoxArray(cba.minimalBox());
    root.cdmap[lev] = amrex::DistributionMapping(root.cba[lev]);
    root.cdmap[lev].define(croot_map);
}

void initializeSpectralWorkspaces(
    const crescent::shared::SolverContext& initialization_context)
{
    for (int lev = 0; lev <= initialization_context.finestLevel(); ++lev) {
        initializeLevelSpectralWorkspace(
            initialization_context.levelContext(lev));
    }
}

void finalizeCoreLevelInitialization(
    const crescent::shared::LevelContext& level_context)
{
    initializeCoreDerivedFields(
        level_context);
    level_context.solutionState().pressure[0].setVal(0.0);
}

void finalizeNestedLevelInitialization(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev > 0);
    auto& solution = level_context.solutionState();
    solution.pressure[lev].setVal(0.0);

    if (crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev)) {
        crescent::nested::setup::initializeNestedL1Level(
            level_context);
    }
}

void initializeLevelState(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    prepareLevelWaveState(
        level_context);
    applyInitializationWaveFilter(
        level_context);
    if (lev == 0) {
        finalizeCoreLevelInitialization(level_context);
    } else {
        finalizeNestedLevelInitialization(level_context);
    }
}

void initializeLevelStates(
    const crescent::shared::SolverContext& initialization_context)
{
    for (int lev = 0; lev <= initialization_context.finestLevel(); ++lev) {
        initializeLevelState(
            initialization_context.levelContext(lev));
    }
}

void printRefinementSummary(
    const crescent::shared::RuntimeControlState& runtime,
    crescent::shared::MeshHierarchyView mesh,
    int finest_level)
{
    // These are solver contracts, not logging diagnostics. Validate them even
    // for quiet runs and for parent-only configurations.
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto level_mesh = mesh.levelMesh(lev);
        const amrex::Box box = level_mesh.grids.minimalBox();
        if (box.length(2) != 1) {
            amrex::Abort(
                "CRESCENT requires exactly one z cell on every level; L" +
                std::to_string(lev) + " has nz=" +
                std::to_string(box.length(2)));
        }
        if (!level_mesh.geom.isPeriodic(0) ||
            !level_mesh.geom.isPeriodic(1)) {
            amrex::Abort(
                "CRESCENT spectral evolution requires periodic geometry in x and y");
        }
    }

    if (finest_level <= 0 || runtime.verbosity < 1) {
        return;
    }

    amrex::Print() << "\n=== AMReX Refinement Ratios ===" << std::endl;
    for (int lev = 0; lev < finest_level; ++lev) {
        amrex::Print() << "  Level " << lev << " -> " << (lev + 1)
                       << ": ref_ratio = " << mesh.ref_ratio[lev] << std::endl;
    }
    for (int lev = 1; lev <= finest_level; ++lev) {
        const auto level_mesh = mesh.levelMesh(lev);
        amrex::Box box = level_mesh.grids.minimalBox();
        amrex::Print() << "  Level " << lev << " box: " << box
                       << " (nz=" << box.length(2) << ")" << std::endl;
    }
    amrex::Print() << "================================\n" << std::endl;
}

} // namespace

void constructSolverState(
    crescent::SolverState& state,
    int max_level,
    const std::string& solver_pp_prefix)
{
    const int nlevs_max = max_level + 1;
    auto& core = state.core;
    core.runtime.resizeLevelState(nlevs_max);
    core.wave.resizeLevelState(nlevs_max);
    for (int lev = 0; lev < nlevs_max; ++lev) {
        core.runtime.nsubsteps[lev] = std::max(1, static_cast<int>(std::pow(2, lev)));
    }

    readParameters(state, max_level, solver_pp_prefix);

    core.solution.resize(nlevs_max);
    if (crescent::nested::hasActiveNestedL1(state.nested.config, max_level)) {
        state.nested.runtime.resize(nlevs_max);
    }
    core.root_layouts.resize(nlevs_max);
    core.spectral_workspace.resize(nlevs_max);
}

void initializeScratchLevelStorage(
    int lev,
    crescent::CoreSolverState& core,
    shared::LevelMeshView level_mesh)
{
    core.solution.defineScratchLevelStorage(
        lev,
        level_mesh.grids,
        level_mesh.dmap,
        core.wave.spectral.perturbation_order);
    core.root_layouts.defineLevel(lev, level_mesh.geom);
}

void initializeSolverState(
    const crescent::shared::SolverContext& context)
{
    printRefinementSummary(
        context.runtime(),
        context.meshHierarchy(),
        context.finestLevel());
    initializeSpectralWorkspaces(context);
    initializeLevelStates(context);
}

void rebuildRestartDerivedState(
    const crescent::shared::SolverContext& context)
{
    auto& core = context.core();
    auto& nested = context.nested();
    const auto mesh = context.meshHierarchy();
    const int finest_level = context.finestLevel();

    printRefinementSummary(
        core.runtime,
        mesh,
        finest_level);
    initializeSpectralWorkspaces(
        context);

    for (int lev = 0; lev <= finest_level; ++lev) {
        initializeLevelSpectralMetrics(
            context.levelContext(lev));
        core.solution.eta[lev].FillBoundary(mesh.geom[lev].periodicity());
        core.solution.phi[lev].FillBoundary(mesh.geom[lev].periodicity());
        core.solution.pressure[lev].FillBoundary(mesh.geom[lev].periodicity());
    }

    initializeCoreDerivedFields(
        context.levelContext(0));

    if (!crescent::nested::hasActiveNestedL1(nested.config, finest_level)) {
        return;
    }

    const int lev = crescent::nested::ActiveNestedL1Level;
    auto l1_context = context.levelContext(lev);
    crescent::nested::setup::initializeDomainsAndStorage(l1_context);
    crescent::nested::setup::zeroDerivedFields(l1_context);

    crescent::nested::snapshots::storeCoreInterpolationSourceAtStepEnd(
        finest_level,
        0,
        core,
        nested,
        mesh.geom);

    crescent::nested::free_surface_evolution::recomputeL1OutputDerivatives(
        l1_context);
    crescent::solver::poststep::computeVelocityOutputs(l1_context);
}

} // namespace crescent::solver::setup

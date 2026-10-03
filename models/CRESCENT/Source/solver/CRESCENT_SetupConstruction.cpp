#include "CRESCENT_SetupLifecycle.H"

#include <algorithm>
#include <cmath>
#include <string>

#include <AMReX_Print.H>

#include "CRESCENT_NestedConfig.H"

namespace {

bool hasConfiguredNestedLevel(
    const crescent::nested::AlgorithmConfig& nested_config,
    int max_level)
{
    return crescent::nested::hasActiveNestedL1(
        nested_config,
        max_level);
}

void resizeSolutionLevelContainers(
    crescent::shared::SolutionState& solution,
    int nlevs_max)
{
    solution.resize(nlevs_max);
}

void resizeNestedLevelContainers(
    crescent::nested::RuntimeState& nested_runtime,
    int nlevs_max)
{
    nested_runtime.resize(nlevs_max);
}

void resizeRootLayoutContainers(
    crescent::shared::RootLayouts& root_layouts,
    int nlevs_max)
{
    root_layouts.resize(nlevs_max);
}

void resizeSpectralWorkspaceContainers(
    crescent::shared::SpectralWorkspace& spectral_workspace,
    int nlevs_max)
{
    spectral_workspace.resize(nlevs_max);
}

void initializeConfiguredLevelStorage(
    const crescent::solver::setup::ConstructionContext& construction_context)
{
    const int nlevs_max = construction_context.levelCount();
    auto core = construction_context.coreAccess();
    auto nested = construction_context.nestedAccess();
    resizeSolutionLevelContainers(core.solution, nlevs_max);

    if (hasConfiguredNestedLevel(
            nested.config,
            construction_context.maxLevel())) {
        resizeNestedLevelContainers(nested.runtime, nlevs_max);
    }

    resizeRootLayoutContainers(core.root_layouts, nlevs_max);
    resizeSpectralWorkspaceContainers(core.spectral_workspace, nlevs_max);
}

void initializeProblemState(
    std::unique_ptr<crescent::ProblemBase>& problem,
    const amrex::Geometry& base_geometry)
{
    problem = crescent::amrex_probinit(
        base_geometry.ProbLo(),
        base_geometry.ProbHi());
}

} // namespace

namespace crescent::solver::setup {

void initializeLevelTimeState(
    const ConstructionContext& construction_context)
{
    auto& runtime = construction_context.runtime();
    auto& wave = construction_context.wave();
    const int nlevs_max = construction_context.levelCount();
    runtime.resizeLevelState(nlevs_max);
    wave.resizeLevelState(nlevs_max);

    for (int lev = 0; lev < nlevs_max; ++lev) {
        runtime.nsubsteps[lev] = std::max(1, static_cast<int>(std::pow(2, lev)));
    }
}

void initializeConfiguredConstructionState(
    const ConstructionContext& construction_context)
{
    auto amr = construction_context.amrAccess();
    initializeConfiguredLevelStorage(construction_context);
    initializeProblemState(amr.prob, construction_context.baseGeometry());
}

void printRefinementSummary(
    const shared::RuntimeControlState& runtime,
    shared::MeshHierarchyView mesh,
    int finest_level)
{
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
        if (box.length(2) != 1) {
            amrex::Abort("ERROR: L" + std::to_string(lev) + " has nz=" +
                         std::to_string(box.length(2)) +
                         " instead of 1! Check amr.ref_ratio_vect in input file.");
        }
    }
    amrex::Print() << "================================\n" << std::endl;
}

} // namespace crescent::solver::setup

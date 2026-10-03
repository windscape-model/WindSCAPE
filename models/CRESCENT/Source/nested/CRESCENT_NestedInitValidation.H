#pragma once

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <string>

#include "CRESCENT_LevelInitializationContext.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_NestedDomains.H"

namespace crescent::nested::validation {

struct InitializationValidationContext
{
    int lev;
    crescent::ownership::CoreSolverAccess core;
    crescent::ownership::NestedSolverAccess nested;
    crescent::shared::MeshHierarchyView mesh;
    crescent::shared::LevelWaveOperatorView level_wave_operator;
    const amrex::BoxArray& root_grids;
    const amrex::DistributionMapping& root_dmap;
    const crescent::nested::domains::DomainStructure& domain;
    int c3_buffer_width;
    const amrex::MultiFab& sigma_blend;
    std::string diagnostic_output_file;

    [[nodiscard]] crescent::shared::LevelMeshView coreLevelMesh() const
    {
        return mesh.levelMesh(0);
    }

    [[nodiscard]] crescent::shared::LevelMeshView levelMesh() const
    {
        return mesh.levelMesh(lev);
    }

    [[nodiscard]] amrex::IntVect refineRatio() const
    {
        return mesh.refinementRatio(lev);
    }
};

inline std::string makeInitializationDiagnosticOutputFile(
    const std::string& plot_file)
{
    constexpr const char* diagnostic_filename =
        "CRESCENT_derivative_diagnostic.dat";
    const std::string::size_type separator = plot_file.find_last_of("/\\");
    if (separator == std::string::npos) {
        return diagnostic_filename;
    }

    return plot_file.substr(0, separator + 1) + diagnostic_filename;
}

inline InitializationValidationContext makeInitializationValidationContext(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    const auto nested = level_context.nestedAccess();
    const auto mesh = level_context.meshHierarchy();
    const auto& coupling_zones = nested.runtime.coupling_zones;
    return {
        lev,
        core,
        nested,
        mesh,
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            mesh.levelMesh(lev).domain()),
        core.root_layouts.ba[lev],
        core.root_layouts.dmap[lev],
        coupling_zones.domains[lev],
        coupling_zones.c3_buffer_width_computed[lev],
        coupling_zones.sigma_blend[lev],
        makeInitializationDiagnosticOutputFile(core.runtime.reporting.plot_file)
    };
}

void runInitializationValidation(
    const InitializationValidationContext& context);

} // namespace crescent::nested::validation

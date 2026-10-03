#ifndef CRESCENT_CONSTRUCTION_CONTEXT_H
#define CRESCENT_CONSTRUCTION_CONTEXT_H

#include <string>
#include <utility>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

#include "CRESCENT_LevelMeshView.H"
#include "ProbBase.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverAccess.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_WaveState.H"

namespace crescent::solver::setup {

class ConstructionContext
{
public:
    ConstructionContext(
        int max_level_in,
        int finest_level_in,
        ownership::CoreSolverAccess core_in,
        ownership::NestedSolverAccess nested_in,
        ownership::AmrIntegrationAccess amr_in,
        const amrex::Geometry& base_geometry_in,
        std::string solver_pp_prefix_in = "crescent")
        : max_level_(max_level_in),
          finest_level_(finest_level_in),
          core_(core_in),
          nested_(nested_in),
          amr_(amr_in),
          base_geometry_(base_geometry_in),
          solver_pp_prefix_(std::move(solver_pp_prefix_in))
    {}

    [[nodiscard]] int maxLevel() const noexcept
    {
        return max_level_;
    }

    [[nodiscard]] int finestLevel() const noexcept
    {
        return finest_level_;
    }

    [[nodiscard]] int levelCount() const noexcept
    {
        return max_level_ + 1;
    }

    [[nodiscard]] ownership::CoreSolverAccess coreAccess() const noexcept
    {
        return core_;
    }

    [[nodiscard]] ownership::NestedSolverAccess nestedAccess() const noexcept
    {
        return nested_;
    }

    [[nodiscard]] ownership::AmrIntegrationAccess amrAccess() const noexcept
    {
        return amr_;
    }

    [[nodiscard]] const amrex::Geometry& baseGeometry() const noexcept
    {
        return base_geometry_;
    }

    [[nodiscard]] shared::RuntimeControlState& runtime() const
    {
        return core_.runtime;
    }

    [[nodiscard]] shared::WaveRuntimeState& wave() const
    {
        return core_.wave;
    }

    [[nodiscard]] const std::string& solverParmParsePrefix() const noexcept
    {
        return solver_pp_prefix_;
    }

private:
    int max_level_;
    int finest_level_;
    ownership::CoreSolverAccess core_;
    ownership::NestedSolverAccess nested_;
    ownership::AmrIntegrationAccess amr_;
    const amrex::Geometry& base_geometry_;
    std::string solver_pp_prefix_;
};

class ScratchLevelContext
{
public:
    ScratchLevelContext(
        int lev_in,
        ownership::CoreSolverAccess core_in,
        ownership::AmrIntegrationAccess amr_in,
        shared::LevelMeshView level_mesh_in)
        : lev_(lev_in),
          core_(core_in),
          amr_(amr_in),
          level_mesh_(level_mesh_in)
    {}

    [[nodiscard]] int level() const noexcept
    {
        return lev_;
    }

    [[nodiscard]] ownership::CoreSolverAccess coreAccess() const noexcept
    {
        return core_;
    }

    [[nodiscard]] ownership::AmrIntegrationAccess amrAccess() const noexcept
    {
        return amr_;
    }

    [[nodiscard]] shared::LevelMeshView levelMesh() const noexcept
    {
        return level_mesh_;
    }

    [[nodiscard]] crescent::ProblemBase& problem() const
    {
        return *amr_.prob;
    }

private:
    int lev_;
    ownership::CoreSolverAccess core_;
    ownership::AmrIntegrationAccess amr_;
    shared::LevelMeshView level_mesh_;
};

inline ConstructionContext makeConstructionContext(
    int max_level,
    int finest_level,
    ownership::CoreSolverAccess core,
    ownership::NestedSolverAccess nested,
    ownership::AmrIntegrationAccess amr,
    const amrex::Geometry& base_geometry,
    std::string solver_pp_prefix = "crescent")
{
    return ConstructionContext(
        max_level,
        finest_level,
        core,
        nested,
        amr,
        base_geometry,
        std::move(solver_pp_prefix));
}

inline ScratchLevelContext makeScratchLevelContext(
    int lev,
    ownership::CoreSolverAccess core,
    ownership::AmrIntegrationAccess amr,
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    return ScratchLevelContext(
        lev,
        core,
        amr,
        {ba, dm, geom});
}

} // namespace crescent::solver::setup

#endif // CRESCENT_CONSTRUCTION_CONTEXT_H

#ifndef CRESCENT_LEVEL_EVOLUTION_CONTEXT_H
#define CRESCENT_LEVEL_EVOLUTION_CONTEXT_H

#include <AMReX_IntVect.H>

#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverAccess.H"
#include "CRESCENT_SpectralWorkspace.H"

namespace crescent::shared {

class LevelEvolutionContext
{
public:
    LevelEvolutionContext(
        int lev_in,
        int finest_level_in,
        ownership::CoreSolverAccess core_in,
        ownership::NestedSolverAccess nested_in,
        MeshHierarchyView mesh_in,
        const amrex::IntVect& core_to_l1_ref_ratio_in)
        : lev_(lev_in),
          finest_level_(finest_level_in),
          core_(core_in),
          nested_(nested_in),
          mesh_(mesh_in),
          core_to_l1_ref_ratio_(core_to_l1_ref_ratio_in)
    {}

    [[nodiscard]] int level() const noexcept
    {
        return lev_;
    }

    [[nodiscard]] int finestLevel() const noexcept
    {
        return finest_level_;
    }

    [[nodiscard]] ownership::CoreSolverAccess coreAccess() const noexcept
    {
        return core_;
    }

    [[nodiscard]] ownership::NestedSolverAccess nestedAccess() const noexcept
    {
        return nested_;
    }

    [[nodiscard]] MeshHierarchyView meshHierarchy() const noexcept
    {
        return mesh_;
    }

    [[nodiscard]] const amrex::IntVect& coreToL1RefinementRatio() const noexcept
    {
        return core_to_l1_ref_ratio_;
    }

    [[nodiscard]] LevelMeshView levelMesh() const
    {
        return mesh_.levelMesh(lev_);
    }

    [[nodiscard]] const amrex::Geometry& levelGeometry() const
    {
        return mesh_.geom[lev_];
    }

    [[nodiscard]] const amrex::Geometry& coreLevelGeometry() const
    {
        return mesh_.geom[0];
    }

    [[nodiscard]] shared::RuntimeControlState& runtime() const
    {
        return core_.runtime;
    }

    [[nodiscard]] shared::SolutionState& solutionState() const
    {
        return core_.solution;
    }

    [[nodiscard]] shared::SpectralWorkspace& spectralWorkspace() const
    {
        return core_.spectral_workspace;
    }

    [[nodiscard]] nested::AlgorithmConfig& nestedConfig() const
    {
        return nested_.config;
    }

    [[nodiscard]] nested::RuntimeState& nestedRuntime() const
    {
        return nested_.runtime;
    }

private:
    int lev_;
    int finest_level_;
    ownership::CoreSolverAccess core_;
    ownership::NestedSolverAccess nested_;
    MeshHierarchyView mesh_;
    amrex::IntVect core_to_l1_ref_ratio_;
};

inline LevelEvolutionContext makeLevelEvolutionContext(
    int lev,
    int finest_level,
    ownership::CoreSolverAccess core,
    ownership::NestedSolverAccess nested,
    MeshHierarchyView mesh,
    const amrex::IntVect& core_to_l1_ref_ratio)
{
    return LevelEvolutionContext(
        lev,
        finest_level,
        core,
        nested,
        mesh,
        core_to_l1_ref_ratio);
}

} // namespace crescent::shared

#endif // CRESCENT_LEVEL_EVOLUTION_CONTEXT_H

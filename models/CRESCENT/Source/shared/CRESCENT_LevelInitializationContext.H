#ifndef CRESCENT_LEVEL_INITIALIZATION_CONTEXT_H
#define CRESCENT_LEVEL_INITIALIZATION_CONTEXT_H

#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverAccess.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_WaveState.H"

namespace crescent::shared {

class LevelInitializationContext
{
public:
    LevelInitializationContext(
        int lev_in,
        int finest_level_in,
        ownership::CoreSolverAccess core_in,
        ownership::NestedSolverAccess nested_in,
        MeshHierarchyView mesh_in)
        : lev_(lev_in),
          finest_level_(finest_level_in),
          core_(core_in),
          nested_(nested_in),
          mesh_(mesh_in)
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

    [[nodiscard]] LevelMeshView levelMesh() const
    {
        return mesh_.levelMesh(lev_);
    }

    [[nodiscard]] LevelMeshView coreLevelMesh() const
    {
        return mesh_.levelMesh(0);
    }

    [[nodiscard]] amrex::IntVect refineRatio() const
    {
        return mesh_.refinementRatio(lev_);
    }

    [[nodiscard]] shared::RuntimeControlState& runtime() const
    {
        return core_.runtime;
    }

    [[nodiscard]] shared::WaveRuntimeState& wave() const
    {
        return core_.wave;
    }

    [[nodiscard]] shared::SolutionState& solutionState() const
    {
        return core_.solution;
    }

    [[nodiscard]] shared::RootLayouts& rootLayouts() const
    {
        return core_.root_layouts;
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
};

inline LevelInitializationContext makeLevelInitializationContext(
    int lev,
    int finest_level,
    ownership::CoreSolverAccess core,
    ownership::NestedSolverAccess nested,
    MeshHierarchyView mesh)
{
    return LevelInitializationContext(
        lev,
        finest_level,
        core,
        nested,
        mesh);
}

} // namespace crescent::shared

#endif // CRESCENT_LEVEL_INITIALIZATION_CONTEXT_H

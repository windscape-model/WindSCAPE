#ifndef CRESCENT_SOLVER_CONTEXT_H
#define CRESCENT_SOLVER_CONTEXT_H

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>

#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_SolverState.H"

namespace crescent::shared {

// The solver state and the mesh hierarchy, seen from one level.
class LevelContext
{
public:
    LevelContext(int lev, int finest_level, SolverState& state, MeshHierarchyView mesh)
        : lev_(lev), finest_level_(finest_level), state_(state), mesh_(mesh)
    {}

    [[nodiscard]] int level() const noexcept { return lev_; }
    [[nodiscard]] int finestLevel() const noexcept { return finest_level_; }

    [[nodiscard]] CoreSolverState& core() const noexcept { return state_.core; }
    [[nodiscard]] NestedSolverState& nested() const noexcept { return state_.nested; }
    [[nodiscard]] RuntimeControlState& runtime() const noexcept { return state_.core.runtime; }
    [[nodiscard]] WaveRuntimeState& wave() const noexcept { return state_.core.wave; }
    [[nodiscard]] SolutionState& solutionState() const noexcept { return state_.core.solution; }
    [[nodiscard]] RootLayouts& rootLayouts() const noexcept { return state_.core.root_layouts; }
    [[nodiscard]] SpectralWorkspace& spectralWorkspace() const noexcept
    {
        return state_.core.spectral_workspace;
    }
    [[nodiscard]] nested::AlgorithmConfig& nestedConfig() const noexcept
    {
        return state_.nested.config;
    }
    [[nodiscard]] nested::RuntimeState& nestedRuntime() const noexcept
    {
        return state_.nested.runtime;
    }

    [[nodiscard]] MeshHierarchyView meshHierarchy() const noexcept { return mesh_; }
    [[nodiscard]] LevelMeshView levelMesh() const { return mesh_.levelMesh(lev_); }
    [[nodiscard]] LevelMeshView coreLevelMesh() const { return mesh_.levelMesh(0); }
    [[nodiscard]] const amrex::Geometry& levelGeometry() const { return mesh_.geom[lev_]; }
    [[nodiscard]] const amrex::Geometry& coreLevelGeometry() const { return mesh_.geom[0]; }
    [[nodiscard]] amrex::IntVect refineRatio() const { return mesh_.refinementRatio(lev_); }

    // FFT, spectral workspace, wavenumbers and filter settings of this level
    // and of the parent level.
    [[nodiscard]] LevelWaveOperatorView waveOperator() const
    {
        return makeLevelWaveOperatorView(lev_, state_.core, mesh_.geom[lev_].data());
    }
    [[nodiscard]] LevelWaveOperatorView coreWaveOperator() const
    {
        return makeLevelWaveOperatorView(0, state_.core, mesh_.geom[0].data());
    }

    // Refinement ratio from the parent (L0) to the child (L1).
    [[nodiscard]] amrex::IntVect coreToL1RefinementRatio() const
    {
        return finest_level_ > 0 ? mesh_.ref_ratio[0] : amrex::IntVect::TheUnitVector();
    }

private:
    int lev_;
    int finest_level_;
    SolverState& state_;
    MeshHierarchyView mesh_;
};

// The solver state and the mesh hierarchy, seen from the whole hierarchy.
class SolverContext
{
public:
    SolverContext(SolverState& state, MeshHierarchyView mesh, int finest_level)
        : state_(state), mesh_(mesh), finest_level_(finest_level)
    {}

    [[nodiscard]] int finestLevel() const noexcept { return finest_level_; }

    [[nodiscard]] CoreSolverState& core() const noexcept { return state_.core; }
    [[nodiscard]] NestedSolverState& nested() const noexcept { return state_.nested; }
    [[nodiscard]] RuntimeControlState& runtime() const noexcept { return state_.core.runtime; }
    [[nodiscard]] SolutionState& solutionState() const noexcept { return state_.core.solution; }
    [[nodiscard]] RootLayouts& rootLayouts() const noexcept { return state_.core.root_layouts; }
    [[nodiscard]] nested::AlgorithmConfig& nestedConfig() const noexcept
    {
        return state_.nested.config;
    }
    [[nodiscard]] nested::RuntimeState& nestedRuntime() const noexcept
    {
        return state_.nested.runtime;
    }

    [[nodiscard]] MeshHierarchyView meshHierarchy() const noexcept { return mesh_; }
    [[nodiscard]] const amrex::Geometry& levelGeometry(int lev) const { return mesh_.geom[lev]; }

    [[nodiscard]] LevelContext levelContext(int lev) const
    {
        return {lev, finest_level_, state_, mesh_};
    }

private:
    SolverState& state_;
    MeshHierarchyView mesh_;
    int finest_level_;
};

} // namespace crescent::shared

#endif // CRESCENT_SOLVER_CONTEXT_H

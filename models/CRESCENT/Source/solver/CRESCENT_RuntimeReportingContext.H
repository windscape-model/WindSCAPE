#ifndef CRESCENT_RUNTIME_REPORTING_CONTEXT_H
#define CRESCENT_RUNTIME_REPORTING_CONTEXT_H

#include <AMReX_IntVect.H>

#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverAccess.H"

namespace crescent::solver {

class RuntimeReportingContext
{
public:
    RuntimeReportingContext(
        ownership::CoreSolverAccess core_in,
        ownership::NestedSolverAccess nested_in,
        shared::MeshHierarchyView mesh_in,
        const amrex::IntVect& core_to_l1_ref_ratio_in,
        int finest_level_in)
        : core_(core_in),
          nested_(nested_in),
          mesh_(mesh_in),
          core_to_l1_ref_ratio_(core_to_l1_ref_ratio_in),
          finest_level_(finest_level_in)
    {}

    [[nodiscard]] shared::RuntimeControlState& runtime() const
    {
        return core_.runtime;
    }

    [[nodiscard]] shared::SolutionState& solutionState() const
    {
        return core_.solution;
    }

    [[nodiscard]] shared::RootLayouts& rootLayouts() const
    {
        return core_.root_layouts;
    }

    [[nodiscard]] nested::AlgorithmConfig& nestedConfig() const
    {
        return nested_.config;
    }

    [[nodiscard]] nested::RuntimeState& nestedRuntime() const
    {
        return nested_.runtime;
    }

    [[nodiscard]] shared::MeshHierarchyView meshHierarchy() const
    {
        return mesh_;
    }

    [[nodiscard]] int finestLevel() const
    {
        return finest_level_;
    }

    [[nodiscard]] const amrex::Geometry& levelGeometry(int lev) const
    {
        return mesh_.geom[lev];
    }

    [[nodiscard]] shared::LevelEvolutionContext levelContext(int lev) const
    {
        return shared::makeLevelEvolutionContext(
            lev,
            finest_level_,
            core_,
            nested_,
            mesh_,
            core_to_l1_ref_ratio_);
    }

private:
    ownership::CoreSolverAccess core_;
    ownership::NestedSolverAccess nested_;
    shared::MeshHierarchyView mesh_;
    amrex::IntVect core_to_l1_ref_ratio_;
    int finest_level_;
};

inline RuntimeReportingContext makeRuntimeReportingContext(
    ownership::CoreSolverAccess core,
    ownership::NestedSolverAccess nested,
    shared::MeshHierarchyView mesh,
    const amrex::IntVect& core_to_l1_ref_ratio,
    int finest_level)
{
    return RuntimeReportingContext(
        core,
        nested,
        mesh,
        core_to_l1_ref_ratio,
        finest_level);
}

} // namespace crescent::solver

#endif // CRESCENT_RUNTIME_REPORTING_CONTEXT_H

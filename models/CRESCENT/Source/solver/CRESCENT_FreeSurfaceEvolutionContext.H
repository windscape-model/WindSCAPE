#ifndef CRESCENT_FREE_SURFACE_EVOLUTION_CONTEXT_H
#define CRESCENT_FREE_SURFACE_EVOLUTION_CONTEXT_H

#include <AMReX_IntVect.H>

#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_RuntimeReportingContext.H"
#include "CRESCENT_SolverAccess.H"

namespace crescent::solver {

class FreeSurfaceEvolutionContext
{
public:
    FreeSurfaceEvolutionContext(
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

    [[nodiscard]] int finestLevel() const
    {
        return finest_level_;
    }

    [[nodiscard]] shared::RuntimeControlState& runtime() const
    {
        return core_.runtime;
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

    [[nodiscard]] RuntimeReportingContext reportingContext() const
    {
        return makeRuntimeReportingContext(
            core_,
            nested_,
            mesh_,
            core_to_l1_ref_ratio_,
            finest_level_);
    }

private:
    ownership::CoreSolverAccess core_;
    ownership::NestedSolverAccess nested_;
    shared::MeshHierarchyView mesh_;
    amrex::IntVect core_to_l1_ref_ratio_;
    int finest_level_;
};

inline FreeSurfaceEvolutionContext makeFreeSurfaceEvolutionContext(
    ownership::CoreSolverAccess core,
    ownership::NestedSolverAccess nested,
    shared::MeshHierarchyView mesh,
    const amrex::IntVect& core_to_l1_ref_ratio,
    int finest_level)
{
    return FreeSurfaceEvolutionContext(
        core,
        nested,
        mesh,
        core_to_l1_ref_ratio,
        finest_level);
}

} // namespace crescent::solver

#endif // CRESCENT_FREE_SURFACE_EVOLUTION_CONTEXT_H

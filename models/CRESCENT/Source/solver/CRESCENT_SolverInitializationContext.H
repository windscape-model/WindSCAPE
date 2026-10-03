#ifndef CRESCENT_SOLVER_INITIALIZATION_CONTEXT_H
#define CRESCENT_SOLVER_INITIALIZATION_CONTEXT_H

#include "CRESCENT_LevelInitializationContext.H"

namespace crescent::solver::setup {

class SolverInitializationContext
{
public:
    SolverInitializationContext(
        ownership::CoreSolverAccess core_in,
        ownership::NestedSolverAccess nested_in,
        shared::MeshHierarchyView mesh_in,
        int finest_level_in)
        : core_(core_in),
          nested_(nested_in),
          mesh_(mesh_in),
          finest_level_(finest_level_in)
    {}

    [[nodiscard]] ownership::CoreSolverAccess coreAccess() const noexcept
    {
        return core_;
    }

    [[nodiscard]] ownership::NestedSolverAccess nestedAccess() const noexcept
    {
        return nested_;
    }

    [[nodiscard]] shared::MeshHierarchyView meshHierarchy() const noexcept
    {
        return mesh_;
    }

    [[nodiscard]] int finestLevel() const noexcept
    {
        return finest_level_;
    }

    [[nodiscard]] shared::RuntimeControlState& runtime() const
    {
        return core_.runtime;
    }

    [[nodiscard]] shared::LevelInitializationContext levelContext(int lev) const
    {
        return shared::makeLevelInitializationContext(
            lev,
            finest_level_,
            core_,
            nested_,
            mesh_);
    }

private:
    ownership::CoreSolverAccess core_;
    ownership::NestedSolverAccess nested_;
    shared::MeshHierarchyView mesh_;
    int finest_level_;
};

inline SolverInitializationContext makeSolverInitializationContext(
    ownership::CoreSolverAccess core,
    ownership::NestedSolverAccess nested,
    shared::MeshHierarchyView mesh,
    int finest_level)
{
    return SolverInitializationContext(
        core,
        nested,
        mesh,
        finest_level);
}

} // namespace crescent::solver::setup

#endif // CRESCENT_SOLVER_INITIALIZATION_CONTEXT_H

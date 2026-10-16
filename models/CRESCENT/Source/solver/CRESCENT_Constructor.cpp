#include "CRESCENT.H"
#include "CRESCENT_Setup.H"
#include "CRESCENT_SolverContext.H"

CRESCENT::CRESCENT()
    : solver_state(std::make_unique<crescent::SolverState>())
{
    crescent::solver::setup::constructSolverState(*solver_state, max_level);
}

CRESCENT::CRESCENT(std::string const& pp_prefix)
    : amrex::AmrCore(pp_prefix),
      solver_state(std::make_unique<crescent::SolverState>()),
      m_solver_pp_prefix(pp_prefix)
{
    crescent::solver::setup::constructSolverState(
        *solver_state, max_level, m_solver_pp_prefix);
}

crescent::SolverState& CRESCENT::solverState() noexcept
{
    return *solver_state;
}

const crescent::SolverState& CRESCENT::solverState() const noexcept
{
    return *solver_state;
}

crescent::shared::MeshHierarchyView CRESCENT::meshHierarchy() const
{
    return {grids, dmap, geom, refRatio()};
}

crescent::shared::SolverContext CRESCENT::solverContext()
{
    return {*solver_state, meshHierarchy(), finest_level};
}

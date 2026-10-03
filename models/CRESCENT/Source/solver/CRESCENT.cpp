#include "CRESCENT.H"

#include <AMReX_IntVect.H>

#include "CRESCENT_FreeSurfaceEvolutionContext.H"
#include "CRESCENT_RunLoop.H"
#include "CRESCENT_SolverOwnership.H"

void CRESCENT::Run() {
    const auto solver = solverOwnership().solverAccess();
    const amrex::IntVect core_to_l1_ref_ratio =
        (finest_level > 0) ? refRatio(0) : amrex::IntVect::TheUnitVector();
    const auto& refinement_ratios = refRatio();
    const auto mesh =
        crescent::shared::makeMeshHierarchyView(grids, dmap, geom, refinement_ratios);
    auto evolution_context = crescent::solver::makeFreeSurfaceEvolutionContext(
        solver.core,
        solver.nested,
        mesh,
        core_to_l1_ref_ratio,
        finest_level);
    crescent::solver::run_loop::execute(evolution_context);
}

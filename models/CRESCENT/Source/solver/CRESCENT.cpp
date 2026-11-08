#include "CRESCENT.H"

#include "CRESCENT_RunLoop.H"
#include "CRESCENT_SolverContext.H"

void CRESCENT::Run() {
    auto context = solverContext();
    crescent::solver::run_loop::execute(context);
}

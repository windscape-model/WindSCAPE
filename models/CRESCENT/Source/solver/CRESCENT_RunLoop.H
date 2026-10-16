#ifndef CRESCENT_RUN_LOOP_H
#define CRESCENT_RUN_LOOP_H

namespace crescent::shared {
class SolverContext;
}

namespace crescent::solver::run_loop {

void execute(crescent::shared::SolverContext& evolution_context);

} // namespace crescent::solver::run_loop

#endif // CRESCENT_RUN_LOOP_H

#ifndef CRESCENT_RUN_LOOP_H
#define CRESCENT_RUN_LOOP_H

namespace crescent::solver {

class FreeSurfaceEvolutionContext;

} // namespace crescent::solver

namespace crescent::solver::run_loop {

void execute(crescent::solver::FreeSurfaceEvolutionContext& evolution_context);

} // namespace crescent::solver::run_loop

#endif // CRESCENT_RUN_LOOP_H

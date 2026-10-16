#ifndef CRESCENT_RUNTIME_OUTPUTS_H
#define CRESCENT_RUNTIME_OUTPUTS_H

namespace crescent::shared {
class SolverContext;
}

namespace crescent::solver::runtime_outputs {

void writePlotFile(const crescent::shared::SolverContext& reporting_context);

void writeTecplotFile(const crescent::shared::SolverContext& reporting_context);

} // namespace crescent::solver::runtime_outputs

#endif // CRESCENT_RUNTIME_OUTPUTS_H

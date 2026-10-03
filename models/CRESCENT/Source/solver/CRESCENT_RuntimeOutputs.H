#ifndef CRESCENT_RUNTIME_OUTPUTS_H
#define CRESCENT_RUNTIME_OUTPUTS_H

namespace crescent::solver {
class RuntimeReportingContext;
}

namespace crescent::solver::runtime_outputs {

void writePlotFile(const crescent::solver::RuntimeReportingContext& reporting_context);

void writeTecplotFile(const crescent::solver::RuntimeReportingContext& reporting_context);

} // namespace crescent::solver::runtime_outputs

#endif // CRESCENT_RUNTIME_OUTPUTS_H

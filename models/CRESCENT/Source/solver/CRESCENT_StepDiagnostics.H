#ifndef CRESCENT_STEP_DIAGNOSTICS_H
#define CRESCENT_STEP_DIAGNOSTICS_H

#include <AMReX_REAL.H>

namespace crescent::shared {
class SolverContext;
}

namespace crescent::solver::runtime_reporting {

void initialize(
    crescent::shared::SolverContext& reporting_context,
    amrex::Real time);

void recordPostStep(
    crescent::shared::SolverContext& reporting_context,
    amrex::Real time);

} // namespace crescent::solver::runtime_reporting

#endif // CRESCENT_STEP_DIAGNOSTICS_H

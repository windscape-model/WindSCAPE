#ifndef CRESCENT_STEP_DIAGNOSTICS_H
#define CRESCENT_STEP_DIAGNOSTICS_H

#include <AMReX_Real.H>

namespace crescent::solver {

class RuntimeReportingContext;

} // namespace crescent::solver

namespace crescent::solver::runtime_reporting {

void initialize(
    crescent::solver::RuntimeReportingContext& reporting_context,
    amrex::Real time);

void recordPostStep(
    crescent::solver::RuntimeReportingContext& reporting_context,
    amrex::Real time);

} // namespace crescent::solver::runtime_reporting

#endif // CRESCENT_STEP_DIAGNOSTICS_H

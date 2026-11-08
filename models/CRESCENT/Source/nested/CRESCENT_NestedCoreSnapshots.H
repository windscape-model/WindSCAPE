#ifndef CRESCENT_NESTED_CORE_SNAPSHOTS_H
#define CRESCENT_NESTED_CORE_SNAPSHOTS_H

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include "CRESCENT_SolverState.H"

namespace crescent::nested::snapshots {

void storeCoreInterpolationSourceAtStepEnd(
    int finest_level,
    int lev,
    crescent::CoreSolverState& core,
    crescent::NestedSolverState& nested,
    const amrex::Vector<amrex::Geometry>& geom);

} // namespace crescent::nested::snapshots

#endif // CRESCENT_NESTED_CORE_SNAPSHOTS_H

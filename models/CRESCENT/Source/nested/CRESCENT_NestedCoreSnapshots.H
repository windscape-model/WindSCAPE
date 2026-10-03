#ifndef CRESCENT_NESTED_CORE_SNAPSHOTS_H
#define CRESCENT_NESTED_CORE_SNAPSHOTS_H

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include "CRESCENT_SolverAccess.H"

namespace crescent::nested::snapshots {

void updateCoreInterpolationSourceAtStepStart(
    int finest_level,
    int lev,
    const crescent::ownership::CoreSolverAccess& core,
    const crescent::ownership::NestedSolverAccess& nested,
    const amrex::Vector<amrex::Geometry>& geom);

void storeCoreInterpolationSourceAtStepEnd(
    int finest_level,
    int lev,
    const crescent::ownership::CoreSolverAccess& core,
    const crescent::ownership::NestedSolverAccess& nested,
    const amrex::Vector<amrex::Geometry>& geom);

} // namespace crescent::nested::snapshots

#endif // CRESCENT_NESTED_CORE_SNAPSHOTS_H

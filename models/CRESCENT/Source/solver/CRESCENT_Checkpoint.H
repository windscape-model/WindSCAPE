#pragma once

#include <functional>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace crescent::ownership {
struct SolverOwnership;
}

namespace crescent::solver {
class RuntimeReportingContext;
}

namespace crescent::solver::checkpoint {

struct CheckpointRestartContext
{
    crescent::ownership::SolverOwnership& ownership;
    int max_level;
    amrex::Vector<amrex::BoxArray>& grids;
    amrex::Vector<amrex::DistributionMapping>& dmap;
    const amrex::Vector<amrex::Geometry>& geom;
    const amrex::Vector<amrex::IntVect>& ref_ratio;
    // AMReX hierarchy setters are protected on CRESCENT; checkpoint I/O owns
    // restart state, while the caller owns hierarchy mutation.
    std::function<void(
        int,
        const amrex::Vector<amrex::BoxArray>&,
        const amrex::Vector<amrex::Real>&)> define_hierarchy;
};

void readCheckpointFile(
    const CheckpointRestartContext& restart_context);

void writeCheckpointFile(
    const crescent::solver::RuntimeReportingContext& reporting_context);

} // namespace crescent::solver::checkpoint

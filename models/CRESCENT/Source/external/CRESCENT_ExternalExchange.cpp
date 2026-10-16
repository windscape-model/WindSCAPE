#include "CRESCENT_ExternalExchange.H"

#include <AMReX_Print.H>

#include "CRESCENT.H"
#include "CRESCENT_MeshHierarchyView.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverState.H"

namespace {

void validateExchangeLevel(int lev, int finest_level)
{
    if (lev < 0 || lev > finest_level) {
        amrex::Abort("CRESCENT external surface exchange level is outside the active hierarchy");
    }
}

void validateExchangeState(
    const crescent::shared::SolutionState& solution,
    const crescent::shared::RuntimeControlState& runtime,
    crescent::shared::MeshHierarchyView mesh,
    int lev)
{
    const int required_size = lev + 1;
    if (static_cast<int>(solution.eta.size()) < required_size ||
        static_cast<int>(solution.phi.size()) < required_size ||
        static_cast<int>(solution.pressure.size()) < required_size ||
        static_cast<int>(solution.u.size()) < required_size ||
        static_cast<int>(solution.v.size()) < required_size ||
        static_cast<int>(solution.phi_z.size()) < required_size ||
        static_cast<int>(runtime.istep.size()) < required_size ||
        static_cast<int>(runtime.t_new.size()) < required_size ||
        static_cast<int>(mesh.geom.size()) < required_size ||
        static_cast<int>(mesh.grids.size()) < required_size ||
        static_cast<int>(mesh.dmap.size()) < required_size) {
        amrex::Abort("CRESCENT external surface exchange state is not allocated for the requested level");
    }
}

crescent::external_exchange::SurfaceExchangeMetadata makeSurfaceExchangeMetadata(
    const crescent::shared::RuntimeControlState& runtime,
    int lev,
    int finest_level)
{
    return {
        lev,
        finest_level,
        runtime.istep[lev],
        runtime.t_new[lev]
    };
}

crescent::external_exchange::FreeSurfaceExchangeView makeFreeSurfaceExchangeView(
    crescent::shared::SolutionState& solution,
    const crescent::shared::RuntimeControlState& runtime,
    crescent::shared::MeshHierarchyView mesh,
    int lev,
    int finest_level)
{
    validateExchangeLevel(lev, finest_level);
    validateExchangeState(solution, runtime, mesh, lev);
    return {
        makeSurfaceExchangeMetadata(runtime, lev, finest_level),
        solution.eta[lev],
        solution.phi[lev],
        solution.pressure[lev],
        solution.u[lev],
        solution.v[lev],
        solution.phi_z[lev],
        mesh.geom[lev],
        mesh.grids[lev],
        mesh.dmap[lev]
    };
}

crescent::external_exchange::ConstFreeSurfaceExchangeView makeConstFreeSurfaceExchangeView(
    const crescent::shared::SolutionState& solution,
    const crescent::shared::RuntimeControlState& runtime,
    crescent::shared::MeshHierarchyView mesh,
    int lev,
    int finest_level)
{
    validateExchangeLevel(lev, finest_level);
    validateExchangeState(solution, runtime, mesh, lev);
    return {
        makeSurfaceExchangeMetadata(runtime, lev, finest_level),
        solution.eta[lev],
        solution.phi[lev],
        solution.pressure[lev],
        solution.u[lev],
        solution.v[lev],
        solution.phi_z[lev],
        mesh.geom[lev],
        mesh.grids[lev],
        mesh.dmap[lev]
    };
}

} // namespace

namespace crescent::external_exchange {

ConstFreeSurfaceExchangeView asConstFreeSurfaceExchange(
    FreeSurfaceExchangeView view)
{
    return {
        view.metadata,
        view.eta,
        view.phi,
        view.pressure,
        view.u,
        view.v,
        view.vertical_velocity,
        view.geometry,
        view.box_array,
        view.distribution_mapping
    };
}

} // namespace crescent::external_exchange

int CRESCENT::finestSurfaceExchangeLevel() const
{
    return finest_level;
}

crescent::external_exchange::FreeSurfaceExchangeView CRESCENT::surfaceExchangeView(
    int lev)
{
    auto& solver = solverState();
    const auto mesh = meshHierarchy();
    return makeFreeSurfaceExchangeView(
        solver.core.solution,
        solver.core.runtime,
        mesh,
        lev,
        finest_level);
}

crescent::external_exchange::ConstFreeSurfaceExchangeView CRESCENT::surfaceExchangeView(
    int lev) const
{
    auto& solver = solverState();
    const auto mesh = meshHierarchy();
    return makeConstFreeSurfaceExchangeView(
        solver.core.solution,
        solver.core.runtime,
        mesh,
        lev,
        finest_level);
}

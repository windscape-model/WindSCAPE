#ifndef CRESCENT_FREE_SURFACE_VIEWS_H
#define CRESCENT_FREE_SURFACE_VIEWS_H

#include <AMReX_MultiFab.H>

namespace crescent::shared {

struct FreeSurfaceStateView
{
    amrex::MultiFab& eta;
    amrex::MultiFab& phi;
};

struct ConstFreeSurfaceStateView
{
    const amrex::MultiFab& eta;
    const amrex::MultiFab& phi;
};

inline ConstFreeSurfaceStateView asConstFreeSurfaceState(
    FreeSurfaceStateView state)
{
    return {state.eta, state.phi};
}

} // namespace crescent::shared

#endif // CRESCENT_FREE_SURFACE_VIEWS_H

#ifndef ERF_WAVE_SURFACE_VEL_BC_H_
#define ERF_WAVE_SURFACE_VEL_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

/**
 * View of an externally supplied (wave-model) sea-surface x-velocity for the orbital
 * moving-wall boundary condition (ERF_USE_TERRAIN_VELOCITY builds).
 *
 * ERF owns the two whole-domain x-face 2D slabs (previous / current coupling instant,
 * broadcast to every rank) and this view of them; an external driver (e.g. the WindSCAPE
 * Coupler) refreshes them each coarse step through ERF::setWaveSurfaceVelocity.
 * ERFPhysBCFunct_u::impose_vertical_xvel_bcs reads the view and imposes the surface
 * velocity time-interpolated to the requested fill time, in place of the analytic
 * monochromatic formula. While !active (standalone ERF), the analytic path is used.
 */
struct WaveSurfaceVelBC
{
    const amrex::FArrayBox* u_old = nullptr;  // x-face slab at t_old
    const amrex::FArrayBox* u_new = nullptr;  // x-face slab at t_new
    amrex::Real t_old = amrex::Real(0.0);
    amrex::Real t_new = amrex::Real(0.0);
    bool active = false;
};

#endif

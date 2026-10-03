#pragma once

#include <AMReX_MultiFab.H>
#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace crescent::dynamics {

/**
 * @brief Advance surface using RK4 time integration
 *
 * @param eta Surface elevation (updated)
 * @param phi Velocity potential (updated)
 * @param eta_saved Saved state from start of RK4 step
 * @param phi_saved Saved state from start of RK4 step
 * @param rhs_eta Right-hand side for eta evolution
 * @param rhs_phi Right-hand side for phi evolution
 * @param dt Time step size
 * @param rk_stage RK4 stage (0-3)
 */
void advanceSurfaceRK4(amrex::MultiFab& eta, amrex::MultiFab& phi,
                       amrex::MultiFab& eta_saved, amrex::MultiFab& phi_saved,
                       const amrex::MultiFab& rhs_eta, const amrex::MultiFab& rhs_phi,
                       const amrex::Real dt, const int rk_stage);

/**
 * @brief Compute right-hand side of HOS evolution equations
 * Computes dη/dt and dφ/dt from current state
 *
 * This routine is used by the core single-domain HOS path. The active nested
 * L1 path computes its RHS in Source/nested/CRESCENT_NestedRK4.cpp after C3
 * repair, derivative blending, and one-way core-to-L1 coupling.
 *
 * @param rhs_eta Output: RHS for eta equation
 * @param rhs_phi Output: RHS for phi equation
 * @param eta Surface elevation
 * @param eta_x ∂η/∂x
 * @param eta_y ∂η/∂y
 * @param phi_x ∂φ/∂x
 * @param phi_y ∂φ/∂y
 * @param phi_z ∂φ/∂z at surface
 * @param pressure Atmospheric pressure
 * @param compute_rhs_metric Temporary workspace for metric term
 * @param compute_rhs_phi_z_sq Temporary workspace for φ_z²
 * @param fft FFT object for dealiasing
 * @param geom Geometry data
 * @param inverse_gravity 1/g
 * @param perturbation_order Order of HOS expansion
 * @param rk_stage RK4 stage number
 */
void computeSurfaceRHS(amrex::MultiFab& rhs_eta, amrex::MultiFab& rhs_phi,
                       const amrex::MultiFab& eta,
                       const amrex::MultiFab& eta_x, const amrex::MultiFab& eta_y,
                       const amrex::MultiFab& phi_x, const amrex::MultiFab& phi_y,
                       const amrex::MultiFab& phi_z,
                       const amrex::MultiFab& pressure,
                       amrex::MultiFab& compute_rhs_metric,
                       amrex::MultiFab& compute_rhs_phi_z_sq,
                       std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
                       const amrex::GeometryData& geom,
                       const amrex::Real inverse_gravity,
                       const int perturbation_order,
                       const int rk_stage,
                       amrex::MultiFab& dealiasing_scratch);

/**
 * @brief Update surface fields using RK4 final combination
 * Combines all RK4 stages: state = saved + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
 *
 * @param eta Surface elevation (updated)
 * @param phi Velocity potential (updated)
 * @param eta_saved Initial state
 * @param phi_saved Initial state
 * @param rhs_eta RHS from all 4 RK4 stages
 * @param rhs_phi RHS from all 4 RK4 stages
 * @param dt Time step size
 */
void updateSurface(amrex::MultiFab& eta, amrex::MultiFab& phi,
                   const amrex::MultiFab& eta_saved, const amrex::MultiFab& phi_saved,
                   const amrex::MultiFab& rhs_eta, const amrex::MultiFab& rhs_phi,
                   const amrex::Real dt);

} // namespace crescent::dynamics

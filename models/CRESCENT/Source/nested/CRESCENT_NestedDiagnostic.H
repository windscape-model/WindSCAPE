#ifndef CRESCENT_NESTED_DIAGNOSTIC_H_
#define CRESCENT_NESTED_DIAGNOSTIC_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <string>

namespace crescent::nested::diagnostic {

/**
 * @brief Write CRESCENT derivative pipeline diagnostic to Tecplot format
 *
 * Outputs all intermediate stages of the CRESCENT derivative computation:
 * - η_interp: surface elevation interpolated from L0
 * - η_periodic/η: surface elevation (with/without C³-repair)
 * - η_x_interp: ∂η/∂x interpolated from L0
 * - η_periodic_x/η_x: ∂η/∂x from FFT on prepared field
 * - φ_interp: interpolated from L0
 * - φ_periodic: after C³-repair
 * - φ_x_interp: ∂φ/∂x interpolated from L0
 * - φ_periodic_x: ∂φ/∂x from FFT on the C³-repaired field
 * - φ_x_blended: ∂φ/∂x after blending (final result)
 * - φ_*_y : same for y-direction (should be ~0 for x-direction wave)
 * - φ_z_L1: ∂φ/∂z computed on L1 using DNO
 * - φ_z_L0_interp: ∂φ/∂z interpolated from L0
 * - φ_z_blended: ∂φ/∂z after blending
 */
void writeDiagnosticTecplot(
    const std::string& filename,
    const amrex::Geometry& geom,
    const amrex::MultiFab& eta_interp,
    const amrex::MultiFab& eta_periodic,
    const amrex::MultiFab& eta_x_interp,
    const amrex::MultiFab& eta_periodic_x,
    const amrex::MultiFab& eta_y_interp,      // ∂η/∂y interpolated from L0
    const amrex::MultiFab& eta_periodic_y,
    const amrex::MultiFab& phi_interp,
    const amrex::MultiFab& phi_periodic,
    const amrex::MultiFab& phi_x_interp,
    const amrex::MultiFab& phi_periodic_x,
    const amrex::MultiFab& phi_x_blended,
    const amrex::MultiFab& phi_y_interp,
    const amrex::MultiFab& phi_periodic_y,
    const amrex::MultiFab& phi_y_blended,
    const amrex::MultiFab& phi_z_L1,
    const amrex::MultiFab& phi_z_L0_interp,
    const amrex::MultiFab& phi_z_blended,
    const amrex::BoxArray& ba_root,
    const amrex::DistributionMapping& dmap_root
);

} // namespace crescent::nested::diagnostic

#endif // CRESCENT_NESTED_DIAGNOSTIC_H_

#ifndef CRESCENT_PHI_EXPANSION_RESEARCH_DIAGNOSTICS_H
#define CRESCENT_PHI_EXPANSION_RESEARCH_DIAGNOSTICS_H

#include <algorithm>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_PhysicsResearchDiagnostics.H"

namespace crescent::physics::research_diagnostics {

inline void logSurfacePotentialSpectrum(
    const amrex::cMultiFab& phi_component_spectral,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    const auto extrema = measureSpectralEnergyExtrema(phi_component_spectral);
    amrex::Print() << "  phi_comp spectral (after FFT): j=0 max=" << extrema.j0_max
                   << ", j≠0 max=" << extrema.other_max
                   << ", ratio="
                   << (extrema.other_max /
                       std::max(extrema.j0_max, static_cast<amrex::Real>(1.0e-30)))
                   << std::endl;
}

inline void logFirstExpansionTransferSpectrum(
    const amrex::cMultiFab& deriv_spectral,
    int kw,
    int k1,
    int verbosity)
{
    if (!(shouldLogDetailedDiagnostics(verbosity) && kw == 1 && k1 == 0)) {
        return;
    }

    const auto extrema = measureSpectralEnergyExtrema(deriv_spectral);
    amrex::Print() << "  dr_c spectral (after multiply, kw=1, k1=0): j=0 max="
                   << extrema.j0_max
                   << ", j≠0 max=" << extrema.other_max
                   << ", ratio="
                   << (extrema.other_max /
                       std::max(extrema.j0_max, static_cast<amrex::Real>(1.0e-30)))
                   << std::endl;
}

inline void logFirstExpansionYUniformity(
    const amrex::MultiFab& expansion_deriv,
    const amrex::MultiFab& eta_pow,
    const amrex::GeometryData& geom,
    int kw,
    int k1,
    int verbosity)
{
    if (!(shouldLogDetailedDiagnostics(verbosity) && kw == 1 && k1 == 0)) {
        return;
    }

    const auto deriv_extrema =
        measureCenterColumnExtrema(expansion_deriv, geom, 0);
    amrex::Print() << "  expansion_deriv (dr, after IFFT, kw=1, k1=0) "
                   << "center-column y-variation (i=" << deriv_extrema.i
                   << ", x=" << deriv_extrema.x << "): min="
                   << deriv_extrema.extrema.min
                   << ", max=" << deriv_extrema.extrema.max
                   << ", Δy="
                   << (deriv_extrema.extrema.max - deriv_extrema.extrema.min)
                   << std::endl;

    const auto eta_extrema =
        measureCenterColumnExtrema(eta_pow, geom, k1);
    amrex::Print() << "  eta_pow[" << k1 << "] (kw=1, k1=0) "
                   << "center-column y-variation (i=" << eta_extrema.i
                   << ", x=" << eta_extrema.x << "): min="
                   << eta_extrema.extrema.min
                   << ", max=" << eta_extrema.extrema.max
                   << ", Δy="
                   << (eta_extrema.extrema.max - eta_extrema.extrema.min)
                   << std::endl;
}

} // namespace crescent::physics::research_diagnostics

#endif // CRESCENT_PHI_EXPANSION_RESEARCH_DIAGNOSTICS_H

#ifndef CRESCENT_PHYSICS_INITIALIZATION_LOGGING_H
#define CRESCENT_PHYSICS_INITIALIZATION_LOGGING_H

#include "CRESCENT_PhysicsResearchDiagnostics.H"
#include "CRESCENT_Spectral.H"

namespace crescent::physics::initialization_logging {

inline void logSurfaceDerivativeSummary(
    const amrex::MultiFab& eta,
    const amrex::MultiFab& eta_x,
    const amrex::MultiFab& eta_y,
    int verbosity)
{
    if (!crescent::physics::research_diagnostics::shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "After computeSurfaceDerivatives, before applyLowPassFilter:"
                   << std::endl;
    amrex::Print() << "  eta min/max: " << eta.min(0) << " / "
                   << eta.max(0) << std::endl;
    amrex::Print() << "  eta_x min/max: " << eta_x.min(0) << " / "
                   << eta_x.max(0) << std::endl;
    amrex::Print() << "  eta_y min/max: " << eta_y.min(0) << " / "
                   << eta_y.max(0) << std::endl;

    amrex::Print() << "After computeSurfaceDerivatives:" << std::endl;
    amrex::Print() << "  eta min/max: " << eta.min(0) << " / "
                   << eta.max(0) << std::endl;
}

inline void logPotentialSpectrum(
    const amrex::MultiFab& phi,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const amrex::GeometryData& geom,
    int verbosity)
{
    if (!crescent::physics::research_diagnostics::shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    auto const& [cba, cdm] = fft->getSpectralDataLayout();
    amrex::cMultiFab phi_spectral(cba, cdm, 1, 0);
    crescent::spectral::forwardFFT(fft, geom, phi, phi_spectral);

    const auto extrema =
        crescent::physics::research_diagnostics::measureSpectralEnergyExtrema(
            phi_spectral);
    amrex::Print() << "  phi spectral energy: j=0 (ky=0) max="
                   << extrema.j0_max
                   << ", j≠0 max=" << extrema.other_max
                   << ", ratio="
                   << (extrema.other_max /
                       std::max(extrema.j0_max, static_cast<amrex::Real>(1.0e-30)))
                   << std::endl;
}

inline void logPhiExpansionState(
    const amrex::MultiFab& phi_expansion,
    const amrex::MultiFab& eta_pow,
    int perturbation_order,
    const amrex::GeometryData& geom,
    int verbosity)
{
    if (!crescent::physics::research_diagnostics::shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "After computeEtaPowers: eta_pow min/max = "
                   << eta_pow.min(0) << " / " << eta_pow.max(0) << std::endl;
    amrex::Print() << "After expandVelocityPotential:" << std::endl;
    amrex::Print() << "  phi_expansion[0] min/max = "
                   << phi_expansion.min(0) << " / " << phi_expansion.max(0)
                   << std::endl;
    const auto extrema =
        crescent::physics::research_diagnostics::measureCenterColumnExtrema(
            phi_expansion,
            geom,
            0);
    amrex::Print() << "  phi_expansion[0] center-column y-variation (i="
                   << extrema.i << ", x=" << extrema.x << "): Δ="
                   << (extrema.extrema.max - extrema.extrema.min) << std::endl;

    if (perturbation_order > 1) {
        amrex::Print() << "  phi_expansion[1] min/max = "
                       << phi_expansion.min(1) << " / "
                       << phi_expansion.max(1) << std::endl;
    }
    if (perturbation_order > 2) {
        amrex::Print() << "  phi_expansion[2] min/max = "
                       << phi_expansion.min(2) << " / "
                       << phi_expansion.max(2) << std::endl;
    }

    amrex::Print() << "  eta_pow[0] max = " << eta_pow.max(0) << std::endl;
    if (perturbation_order > 1) {
        amrex::Print() << "  eta_pow[1] max = " << eta_pow.max(1) << std::endl;
    }
}

inline void logVerticalVelocityState(
    const amrex::MultiFab& phi_z,
    const amrex::GeometryData& geom,
    int verbosity)
{
    if (verbosity >= 1) {
        amrex::Print()
            << "After computeVerticalVelocity (initialization): phi_z min/max = "
            << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
    }

    if (!crescent::physics::research_diagnostics::shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    const auto extrema =
        crescent::physics::research_diagnostics::measureCenterColumnExtrema(
            phi_z,
            geom,
            0);
    amrex::Print() << "  L0 phi_z center-column y-variation at initialization (i="
                   << extrema.i << ", x=" << extrema.x << "): Δ="
                   << (extrema.extrema.max - extrema.extrema.min) << std::endl;
}

} // namespace crescent::physics::initialization_logging

#endif // CRESCENT_PHYSICS_INITIALIZATION_LOGGING_H

#ifndef CRESCENT_LEVEL_WAVE_OPERATOR_VIEW_H
#define CRESCENT_LEVEL_WAVE_OPERATOR_VIEW_H

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

#include "CRESCENT_SolverAccess.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_WaveState.H"

namespace crescent::shared {

struct LevelWaveOperatorView
{
    int lev;
    int perturbation_order;
    amrex::GeometryData geomdata;
    SpectralWorkspace& spectral_workspace;
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft;
    LevelWaveSpectralSettings spectral_settings;

    [[nodiscard]] const WavenumberSpacing& wavenumbers() const
    {
        return spectral_settings.wavenumbers;
    }

    [[nodiscard]] const SpectralFilterSettings& filterSettings() const
    {
        return spectral_settings.filter;
    }
};

inline LevelWaveOperatorView makeLevelWaveOperatorView(
    int lev,
    const WaveRuntimeState& wave,
    SpectralWorkspace& spectral_workspace,
    amrex::GeometryData geomdata)
{
    return {
        lev,
        wave.spectral.perturbation_order,
        geomdata,
        spectral_workspace,
        spectral_workspace.fft[lev],
        wave.levelSpectralSettings(lev)
    };
}

inline LevelWaveOperatorView makeLevelWaveOperatorView(
    int lev,
    const crescent::ownership::CoreSolverAccess& core,
    amrex::GeometryData geomdata)
{
    return makeLevelWaveOperatorView(
        lev,
        core.wave,
        core.spectral_workspace,
        geomdata);
}

} // namespace crescent::shared

#endif // CRESCENT_LEVEL_WAVE_OPERATOR_VIEW_H

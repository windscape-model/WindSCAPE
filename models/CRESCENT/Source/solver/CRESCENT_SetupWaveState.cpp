#include "CRESCENT_SetupWaveState.H"

#include "CRESCENT_CoreWaveInitialization.H"

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_LevelInitLogging.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_SpectralGridMetrics.H"
#include "CRESCENT_Spectral.H"

namespace {

void computeLevelWavenumberSpacing(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    auto& wave = level_context.wave();
    const int verbosity = level_context.runtime().verbosity;
    const auto level_mesh = level_context.levelMesh();
    const auto spectral_grid =
        crescent::shared::makeSpectralGridMetrics(
            level_mesh.grids.minimalBox(),
            level_mesh.geom);
    auto& wavenumbers = wave.levelWavenumbers(lev);
    wavenumbers = spectral_grid.wavenumberSpacing();

    if (verbosity >= 1) {
        amrex::Print() << "Computed wavenumber spacing: dk_x[" << lev << "] = "
                       << wavenumbers.dk_x
                       << ", dk_y[" << lev << "] = "
                       << wavenumbers.dk_y
                       << " (Lx=" << spectral_grid.Lx
                       << ", Ly=" << spectral_grid.Ly << ")" << std::endl;
    }
}

void computeLevelWavenumberPowers(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    const auto& wave = level_context.wave();
    const auto& wavenumbers = wave.levelWavenumbers(lev);
    crescent::spectral::computeWavenumberPowers(
        core.solution.wavenumber_pow[lev],
        wavenumbers.dk_x,
        wavenumbers.dk_y);
    if (level_context.runtime().verbosity >= 2) {
        amrex::Print() << "After calculate_wavenum:" << std::endl;
    }
    crescent::solver::logging::logWavenumberPowerSamples(
        core.solution, core.runtime, core.wave, lev);
}

} // namespace

namespace crescent::solver::setup::wave_state {

void initializeLevelSpectralMetrics(
    const crescent::shared::LevelInitializationContext& level_context)
{
    computeLevelWavenumberSpacing(level_context);
    computeLevelWavenumberPowers(level_context);
}

void prepareLevelWaveState(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    computeLevelWavenumberSpacing(level_context);

    if (lev == 0) {
        crescent::init::initializeCoreWaveField(
            core,
            level_context.levelMesh(),
            level_context.finestLevel());
    } else {
        AMREX_ALWAYS_ASSERT(lev > 0);
        if (level_context.runtime().verbosity >= 1) {
            amrex::Print() << "\n=== Skipping wave initialization for L" << lev
                           << " (will be interpolated from parent) ===" << std::endl;
        }
    }

    computeLevelWavenumberPowers(level_context);
}

void applyInitializationWaveFilter(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            level_context.levelMesh().domain());
    auto free_surface_state = core.solution.freeSurfaceStateAt(lev);
    crescent::spectral::applyLowPassFilter(
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        free_surface_state,
        level_wave_operator.spectral_workspace.eta_spectral[lev],
        level_wave_operator.spectral_settings,
        level_wave_operator.spectral_workspace.peak_wavenumber_gather_scratch[lev],
        level_wave_operator.spectral_workspace.peak_wavenumber_reduction_scratch[lev]);
}

} // namespace crescent::solver::setup::wave_state

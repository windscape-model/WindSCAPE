#ifndef CRESCENT_WAVE_STATE_H
#define CRESCENT_WAVE_STATE_H

#include <string>

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace crescent::shared {

struct WavenumberSpacing
{
    amrex::Real dk_x = 0.0;
    amrex::Real dk_y = 0.0;
};

struct SpectralFilterSettings
{
    amrex::Real peak_wavenumber_physical = -1.0;
    amrex::Real beta1_filter = 8.0;
    amrex::Real beta2_filter = 12.0;
};

struct LevelWaveSpectralSettings
{
    WavenumberSpacing wavenumbers;
    SpectralFilterSettings filter;
};

struct WaveInitializationState
{
    amrex::Real wind_speed_10m = 0.0;
    amrex::Real fetch = 0.0;
    amrex::Real alpha = 0.0;
    amrex::Real peak_frequency = 1.0;
    amrex::Real peak_wavelength = 0.0;
    amrex::Real wave_steepness = 0.1;
    amrex::Real swell_wavenumber = 1.0;
    amrex::Real initial_angle = 0.0;
    amrex::Real angle_spread = 1.0;

    std::string wave_type = "jonswap";  // Options: "stokes", "jonswap", "monochromatic", "focused_group"
    amrex::Long jonswap_random_seed = 12345; // Fixed JONSWAP phase seed for reproducible initialization (-1 = use AMReX default seeding)

    int directional_exponent = 2;       // Directional spreading exponent for JONSWAP (typical value: 2)
    int use_directional_spectrum = 1;   // Use directional spreading (1=yes, 0=no)
    int wave_mode_x = 5;                // Wave mode index for Stokes wave - sets total wavenumber magnitude (default 5)
    amrex::Real wave_angle = 0.0;       // Stokes wave propagation angle in radians (0.0 = +x direction, π/2 = +y direction)
};

struct WaveSpectralState
{
    int perturbation_order = 3; // maximum perturbation order for HOS expansion
    amrex::Vector<WavenumberSpacing> wavenumber_spacing; // spectral spacing per level
    amrex::Real beta1_filter = 8.0, beta2_filter = 12.0;
    amrex::Real peak_wavenumber_override = -1.0;  // User-specified peak wavenumber for filter (-1.0 = compute automatically)
    amrex::Real peak_wavenumber_physical = -1.0;  // Physical peak wavenumber (rad/m) computed from wave initialization
};

struct WaveRuntimeState
{
    WaveInitializationState initialization;
    WaveSpectralState spectral;

    void resizeLevelState(int nlevs_max)
    {
        spectral.wavenumber_spacing.resize(nlevs_max);
    }

    [[nodiscard]] WavenumberSpacing& levelWavenumbers(int lev)
    {
        return spectral.wavenumber_spacing[lev];
    }

    [[nodiscard]] const WavenumberSpacing& levelWavenumbers(int lev) const
    {
        return spectral.wavenumber_spacing[lev];
    }

    [[nodiscard]] SpectralFilterSettings spectralFilterSettings() const
    {
        return {
            spectral.peak_wavenumber_physical,
            spectral.beta1_filter,
            spectral.beta2_filter
        };
    }

    [[nodiscard]] LevelWaveSpectralSettings levelSpectralSettings(int lev) const
    {
        return {levelWavenumbers(lev), spectralFilterSettings()};
    }
};

} // namespace crescent::shared

#endif // CRESCENT_WAVE_STATE_H

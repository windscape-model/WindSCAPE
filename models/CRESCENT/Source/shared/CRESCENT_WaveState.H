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
    amrex::Real beta1_filter = 9.0;
    amrex::Real beta2_filter = 12.0;
    // Child (L1) spectral roll-off exp(-alpha (|k|/k_N)^p) applied by the
    // combined L1 filter; the historical values (33000, 20) roll off from
    // ~0.5 k_N (weight 0.97) to ~0.6 k_N (weight 0.30) and are far more
    // selective than the Hou-Li (2007) filter (36, 36).
    amrex::Real l1_rolloff_alpha = 33000.0;
    int l1_rolloff_power = 20;
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
    amrex::Real wave_steepness = 0.1;     // Stokes/crescent amplitude control; not used by JONSWAP
    amrex::Real initial_angle = 0.0;
    amrex::Real angle_spread = 1.0;

    std::string wave_type = "jonswap";  // Options: "stokes", "crescent", "jonswap"
    amrex::Long jonswap_random_seed = 12345; // Fixed JONSWAP phase seed for reproducible initialization (-1 = use AMReX default seeding)

    int directional_exponent = 2;       // Cosine-power exponent for JONSWAP directional spreading
    int use_directional_spectrum = 1;   // 1=directional spreading, 0=long-crested +x JONSWAP
    amrex::Real jonswap_amplitude_scale = 1.0; // Multiplier on the spectrum-derived JONSWAP Fourier amplitudes
    int wave_mode_x = 5;                // Absolute Stokes x-mode away from the y axis (default 5)
    amrex::Real wave_angle = 0.0;       // Stokes wave propagation angle in radians (0.0 = +x direction, π/2 = +y direction)

    // Crescent (3D horseshoe) initialization: a 5th-order Stokes carrier along
    // +x plus an oblique free-wave perturbation (the class-II instability seed).
    // Ported from the reference Fortran HOS crescent_init; used only when
    // wave_type == "crescent". Carrier k and steepness come from wave_mode_x /
    // wave_steepness; the ratios are the perturbation wavenumbers in carrier units.
    amrex::Real crescent_eps = 0.16;       // perturbation amplitude relative to carrier amplitude a
    amrex::Real crescent_kx_ratio = 1.5;   // perturbation kx / carrier k (streamwise 3/2 subharmonic)
    amrex::Real crescent_ky_ratio = 1.23;  // perturbation ky / carrier k (transverse oblique mode)
    int crescent_stokes_order = 5;         // Stokes order of the analytic carrier (1..5)
};

struct WaveSpectralState
{
    int perturbation_order = 3; // maximum perturbation order for HOS expansion
    amrex::Vector<WavenumberSpacing> wavenumber_spacing; // spectral spacing per level
    amrex::Real beta1_filter = 9.0, beta2_filter = 12.0; // low-pass (5): default beta1 = 9 = the manuscript baseline since 2026-08-16 (was 8)
    amrex::Real beta1_filter_l1 = -1.0; // L1 low-pass width override (-1 = use beta1_filter); lets each level damp its own cutoff-pileup band
    amrex::Real l1_rolloff_alpha = 33000.0; // crescent.l1_rolloff_alpha: L1 spectral roll-off strength exp(-alpha (|k|/k_N)^p) (historical 33000; Hou-Li 2007 uses 36)
    int l1_rolloff_power = 20;              // crescent.l1_rolloff_power: L1 spectral roll-off exponent p (historical 20; Hou-Li 2007 uses 36)
    amrex::Real adjustment_periods = 0.0; // Dommermuth (2000) nonlinear ramp time Ta in peak periods (0 = off); required to linearly initialize steep seas
    amrex::Real peak_wavenumber_override = -1.0;  // User-specified peak wavenumber for filter (-1.0 = use the initializer's peak)
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
        const amrex::Real filter_peak =
            spectral.peak_wavenumber_override > 0.0
            ? spectral.peak_wavenumber_override
            : spectral.peak_wavenumber_physical;
        return {
            filter_peak,
            spectral.beta1_filter,
            spectral.beta2_filter,
            spectral.l1_rolloff_alpha,
            spectral.l1_rolloff_power
        };
    }

    [[nodiscard]] LevelWaveSpectralSettings levelSpectralSettings(int lev) const
    {
        auto filter = spectralFilterSettings();
        if (lev >= 1 && spectral.beta1_filter_l1 > 0.0) {
            filter.beta1_filter = spectral.beta1_filter_l1;
        }
        return {levelWavenumbers(lev), filter};
    }
};

} // namespace crescent::shared

#endif // CRESCENT_WAVE_STATE_H

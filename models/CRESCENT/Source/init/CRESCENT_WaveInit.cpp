#include "CRESCENT_CoreWaveInitialization.H"
#include "CRESCENT_Init.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_WaveState.H"
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>
#include <string>

using namespace amrex;
using namespace crescent::constants;

namespace {

struct StokesWavevector
{
    int mode_x = 0;
    int mode_y = 0;
    Real kx = 0.0;
    Real ky = 0.0;
    Real magnitude = 0.0;
    Real actual_angle = 0.0;
};

StokesWavevector computeStokesWavevector(
    const Real dk_x,
    const Real dk_y,
    const int wave_mode_x,
    const Real wave_angle)
{
    StokesWavevector wavevector;

    if (std::abs(wave_angle) < 1.0e-10) {
        wavevector.mode_x = wave_mode_x;
    } else if (std::abs(wave_angle - PI / 2.0) < 1.0e-10) {
        wavevector.mode_y = wave_mode_x;
    } else {
        wavevector.mode_x = wave_mode_x;
        const Real target_ratio =
            std::cos(wave_angle) / std::sin(wave_angle) * (dk_y / dk_x);
        wavevector.mode_y = std::max(1, int(std::round(wavevector.mode_x / target_ratio)));
    }

    wavevector.kx = wavevector.mode_x * dk_x;
    wavevector.ky = wavevector.mode_y * dk_y;
    wavevector.magnitude = std::sqrt(
        wavevector.kx * wavevector.kx + wavevector.ky * wavevector.ky);
    wavevector.actual_angle = std::atan2(wavevector.ky, wavevector.kx);
    return wavevector;
}

void initializeCoreStokesWave(
    crescent::ownership::CoreSolverAccess core,
    crescent::shared::LevelMeshView level_mesh)
{
    const auto& wavenumbers = core.wave.levelWavenumbers(0);
    const auto& initialization = core.wave.initialization;
    auto& spectral = core.wave.spectral;
    auto free_surface_state = core.solution.freeSurfaceStateAt(0);

    crescent::init::initializeStokesWave(
        free_surface_state,
        core.root_layouts.ba[0],
        core.root_layouts.dmap[0],
        level_mesh.domain(),
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        initialization.wave_steepness,
        crescent::constants::INV_GRAVITY,
        initialization.wave_mode_x,
        initialization.wave_angle);

    const auto wavevector = computeStokesWavevector(
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        initialization.wave_mode_x,
        initialization.wave_angle);
    spectral.peak_wavenumber_physical = wavevector.magnitude;

    amrex::Print() << "Computed physical peak wavenumber: "
                   << spectral.peak_wavenumber_physical
                   << " rad/m (kx=" << wavevector.kx
                   << ", ky=" << wavevector.ky << ")"
                   << std::endl;
}

void configureCoreJonswapTimeStep(
    crescent::shared::RuntimeControlState& runtime,
    const crescent::shared::WaveRuntimeState& wave,
    int finest_level)
{
    const auto& initialization = wave.initialization;
    const Real T_peak = crescent::constants::PI2 / initialization.peak_frequency;
    const Real dt_auto = T_peak / 20.0;

    if (runtime.dt_override > 0.0) {
        runtime.dt[0] = runtime.dt_override;
        amrex::Print() << "\n=== JONSWAP dt (USER OVERRIDE) ===" << std::endl;
        amrex::Print() << "  peak_frequency (rad/s) = "
                       << initialization.peak_frequency << std::endl;
        amrex::Print() << "  T_peak (s) = " << T_peak << std::endl;
        amrex::Print() << "  dt (automatic) = " << dt_auto << " s" << std::endl;
        amrex::Print() << "  dt (user override) = " << runtime.dt[0] << " s" << std::endl;
    } else {
        runtime.dt[0] = dt_auto;
        amrex::Print() << "\n=== JONSWAP Automatic dt Calculation ===" << std::endl;
        amrex::Print() << "  peak_frequency (rad/s) = "
                       << initialization.peak_frequency << std::endl;
        amrex::Print() << "  T_peak (s) = " << T_peak << std::endl;
        amrex::Print() << "  dt = T_peak/20 = " << runtime.dt[0] << " s" << std::endl;
    }
    amrex::Print() << "========================================\n" << std::endl;

    for (int lev = 1; lev <= finest_level; ++lev) {
        runtime.dt[lev] = runtime.dt[0];
        amrex::Print() << "  dt[" << lev << "] = " << runtime.dt[lev]
                       << " s (same as L0)" << std::endl;
    }
}

void initializeCoreJonswapWaves(
    crescent::ownership::CoreSolverAccess core,
    crescent::shared::LevelMeshView level_mesh)
{
    const auto& wavenumbers = core.wave.levelWavenumbers(0);
    auto& initialization = core.wave.initialization;
    auto& spectral = core.wave.spectral;
    auto free_surface_state = core.solution.freeSurfaceStateAt(0);

    if (initialization.jonswap_random_seed >= 0) {
        const amrex::ULong cpu_seed =
            static_cast<amrex::ULong>(initialization.jonswap_random_seed);
        const amrex::ULong gpu_seed = cpu_seed + 987654321ULL;
        amrex::ResetRandomSeed(cpu_seed, gpu_seed);
        amrex::Print() << "Using fixed JONSWAP random seed: cpu=" << cpu_seed
                       << ", gpu=" << gpu_seed << std::endl;
    } else {
        amrex::Print() << "Using AMReX default JONSWAP random seeding" << std::endl;
    }

    crescent::init::generateJonswapWaves(
        free_surface_state,
        core.spectral_workspace.fft[0],
        level_mesh.domain(),
        core.root_layouts.cba[0],
        core.root_layouts.cdmap[0],
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        initialization.wind_speed_10m,
        initialization.fetch,
        crescent::constants::INV_GRAVITY,
        initialization.initial_angle,
        initialization.angle_spread,
        initialization.directional_exponent,
        initialization.alpha,
        initialization.peak_frequency,
        initialization.peak_wavelength,
        core.runtime.verbosity);

    spectral.peak_wavenumber_physical =
        initialization.peak_frequency * initialization.peak_frequency
        * crescent::constants::INV_GRAVITY;

    amrex::Print() << "Computed physical peak wavenumber: "
                   << spectral.peak_wavenumber_physical
                   << " rad/m (from peak_frequency=" << initialization.peak_frequency
                   << " rad/s)" << std::endl;
}

} // namespace

namespace crescent::init {

void initializeStokesWave(MultiFab& eta,
                          MultiFab& phi,
                          const BoxArray& ba_root,
                          const DistributionMapping& dmap_root,
                          const GeometryData& geom,
                          const Real dk_x,
                          const Real dk_y,
                          const Real wave_steepness,
                          const Real inverse_gravity,
                          const int wave_mode_x,
                          const Real wave_angle)
{
    amrex::Print() << "\n=== STOKES_2D INITIALIZATION (ROTATED) ===" << std::endl;
    amrex::Print() << "Parameters: wave_mode_x=" << wave_mode_x << ", wave_angle=" << wave_angle
                   << " rad (" << (wave_angle * 180.0 / PI) << " deg)" << std::endl;
    amrex::Print() << "  dk_x=" << dk_x << ", dk_y=" << dk_y
                   << ", wave_steepness=" << wave_steepness
                   << ", inverse_gravity=" << inverse_gravity << std::endl;

    MultiFab eta_o(ba_root, dmap_root, 1, 0);
    MultiFab phi_o(ba_root, dmap_root, 1, 0);

    eta_o.setVal(0.0);
    phi_o.setVal(0.0);

    if (ParallelDescriptor::IOProcessor()) {
        Real theta = 0.0;

        const auto wavevector =
            computeStokesWavevector(dk_x, dk_y, wave_mode_x, wave_angle);
        const int mode_x = wavevector.mode_x;
        const int mode_y = wavevector.mode_y;
        const Real kx = wavevector.kx;
        const Real ky = wavevector.ky;
        const Real k_mag = wavevector.magnitude;
        const Real actual_angle = wavevector.actual_angle;
        Real omega = std::sqrt(k_mag / inverse_gravity);
        Real amplitude = wave_steepness / k_mag;

        amrex::Print() << "Calculated wave properties:" << std::endl;
        amrex::Print() << "  Mode numbers: (" << mode_x << ", " << mode_y << ")" << std::endl;
        amrex::Print() << "  kx = " << kx << " rad/m (mode_x * dk_x)" << std::endl;
        amrex::Print() << "  ky = " << ky << " rad/m (mode_y * dk_y)" << std::endl;
        amrex::Print() << "  k_mag = " << k_mag << " rad/m" << std::endl;
        amrex::Print() << "  Requested angle = " << (wave_angle * 180.0 / PI) << " deg" << std::endl;
        amrex::Print() << "  Actual angle = " << (actual_angle * 180.0 / PI) << " deg" << std::endl;
        amrex::Print() << "  omega = " << omega << " rad/s" << std::endl;
        amrex::Print() << "  amplitude = " << amplitude << " m" << std::endl;
        amrex::Print() << "Domain: " << geom.Domain() << std::endl;
        amrex::Print() << "dx=" << geom.CellSize()[0] << ", dy=" << geom.CellSize()[1] << std::endl;

        for (MFIter mfi(eta_o, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto eta_o_arr = eta_o.array(mfi);
            auto phi_o_arr = phi_o.array(mfi);

            const auto& dx = geom.CellSize();
            const auto& domain = geom.Domain();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Real x = (i - domain.smallEnd(0)) * dx[0];
                Real y = (j - domain.smallEnd(1)) * dx[1];
                Real phase = kx * x + ky * y + theta;

                eta_o_arr(i,j,k) = amplitude * std::cos(phase);
                phi_o_arr(i,j,k) = (amplitude / inverse_gravity / omega) * std::sin(phase);
            });
        }
    }

    ParallelCopy(eta, eta_o, 0, 0, 1);
    ParallelCopy(phi, phi_o, 0, 0, 1);

    amrex::Print() << "After one_to_all:" << std::endl;
    amrex::Print() << "eta min/max: " << eta.min(0) << " / " << eta.max(0) << std::endl;
    amrex::Print() << "phi min/max: " << phi.min(0) << " / " << phi.max(0) << std::endl;
    amrex::Print() << "=== END STOKES_2D ===" << std::endl << std::endl;
}

void initializeCoreWaveField(
    crescent::ownership::CoreSolverAccess core,
    crescent::shared::LevelMeshView level_mesh,
    int finest_level)
{
    const auto& initialization = core.wave.initialization;
    amrex::Print() << "\n=== Wave Initialization: "
                   << initialization.wave_type << " ===" << std::endl;

    if (initialization.wave_type == "stokes") {
        initializeCoreStokesWave(core, level_mesh);
    } else if (initialization.wave_type == "jonswap") {
        initializeCoreJonswapWaves(core, level_mesh);
        configureCoreJonswapTimeStep(core.runtime, core.wave, finest_level);
    } else {
        amrex::Abort("Unknown wave_type: " + initialization.wave_type +
                     ". Valid options: 'stokes', 'jonswap'");
    }

    const auto free_surface_state =
        crescent::shared::asConstFreeSurfaceState(core.solution.freeSurfaceStateAt(0));

    if (core.runtime.verbosity >= 2) {
        amrex::Print() << "Before copy to eta_saved:" << std::endl;
        amrex::Print() << "  eta min/max: " << free_surface_state.eta.min(0)
                       << " / " << free_surface_state.eta.max(0) << std::endl;
        amrex::Print() << "  phi min/max: " << free_surface_state.phi.min(0)
                       << " / " << free_surface_state.phi.max(0) << std::endl;
    }

    core.solution.saveSurfaceState(0);

    if (core.runtime.verbosity >= 2) {
        amrex::Print() << "After copy to eta_saved:" << std::endl;
        amrex::Print() << "  eta min/max: " << free_surface_state.eta.min(0)
                       << " / " << free_surface_state.eta.max(0) << std::endl;
        amrex::Print() << "  eta_saved min/max: " << core.solution.eta_saved[0].min(0)
                       << " / " << core.solution.eta_saved[0].max(0) << std::endl;
    }
}

} // namespace crescent::init

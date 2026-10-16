#include "CRESCENT_CoreWaveInitialization.H"
#include "CRESCENT_Init.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_StokesWavevector.H"
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

void initializeCoreStokesWave(
    crescent::CoreSolverState& core,
    crescent::shared::LevelMeshView level_mesh)
{
    const auto& wavenumbers = core.wave.levelWavenumbers(0);
    const auto& initialization = core.wave.initialization;
    auto& spectral = core.wave.spectral;

    crescent::init::initializeStokesWave(
        core.solution.eta[0],
        core.solution.phi[0],
        core.root_layouts.ba[0],
        core.root_layouts.dmap[0],
        level_mesh.domain(),
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        initialization.wave_steepness,
        crescent::constants::INV_GRAVITY,
        initialization.wave_mode_x,
        initialization.wave_angle);

    const auto wavevector = crescent::init::stokes::computeWavevector(
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

void initializeCoreCrescentWave(
    crescent::CoreSolverState& core,
    crescent::shared::LevelMeshView level_mesh)
{
    const auto& wavenumbers = core.wave.levelWavenumbers(0);
    const auto& initialization = core.wave.initialization;
    auto& spectral = core.wave.spectral;

    crescent::init::initializeCrescentWave(
        core.solution.eta[0],
        core.solution.phi[0],
        core.root_layouts.ba[0],
        core.root_layouts.dmap[0],
        level_mesh.domain(),
        wavenumbers.dk_x,
        wavenumbers.dk_y,
        initialization.wave_steepness,
        crescent::constants::INV_GRAVITY,
        initialization.wave_mode_x,
        initialization.crescent_eps,
        initialization.crescent_kx_ratio,
        initialization.crescent_ky_ratio,
        initialization.crescent_stokes_order);

    // Filter peak is the carrier magnitude (carrier along +x).
    spectral.peak_wavenumber_physical =
        initialization.wave_mode_x * wavenumbers.dk_x;

    amrex::Print() << "Computed physical peak wavenumber: "
                   << spectral.peak_wavenumber_physical
                   << " rad/m (carrier k)" << std::endl;
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
    crescent::CoreSolverState& core,
    crescent::shared::LevelMeshView level_mesh)
{
    const auto& wavenumbers = core.wave.levelWavenumbers(0);
    auto& initialization = core.wave.initialization;
    auto& spectral = core.wave.spectral;

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
        core.solution.eta[0],
        core.solution.phi[0],
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
        initialization.use_directional_spectrum,
        initialization.jonswap_amplitude_scale,
        initialization.alpha,
        initialization.peak_frequency,
        initialization.peak_wavelength);

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
            crescent::init::stokes::computeWavevector(
                dk_x,
                dk_y,
                wave_mode_x,
                wave_angle);
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

            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx{
                geom.CellSize()[0], geom.CellSize()[1], geom.CellSize()[2]};
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

    amrex::Print() << "=== END STOKES_2D ===" << std::endl << std::endl;
}

void initializeCrescentWave(MultiFab& eta,
                            MultiFab& phi,
                            const BoxArray& ba_root,
                            const DistributionMapping& dmap_root,
                            const GeometryData& geom,
                            const Real dk_x,
                            const Real dk_y,
                            const Real wave_steepness,
                            const Real inverse_gravity,
                            const int wave_mode_x,
                            const Real crescent_eps,
                            const Real crescent_kx_ratio,
                            const Real crescent_ky_ratio,
                            const int crescent_stokes_order)
{
    amrex::Print() << "\n=== CRESCENT (3D horseshoe) INITIALIZATION ===" << std::endl;

    MultiFab eta_o(ba_root, dmap_root, 1, 0);
    MultiFab phi_o(ba_root, dmap_root, 1, 0);

    eta_o.setVal(0.0);
    phi_o.setVal(0.0);

    if (ParallelDescriptor::IOProcessor()) {
        // Carrier: x-propagating (the crescent recipe assumes a 1D carrier).
        const Real kc = wave_mode_x * dk_x;      // carrier wavenumber (|k|)
        const Real ak = wave_steepness;          // carrier steepness = a*kc
        const Real a = ak / kc;                  // carrier amplitude
        const Real g = 1.0 / inverse_gravity;
        const Real sqrt_g_over_k3 = std::sqrt(g / (kc * kc * kc));
        const int order = crescent_stokes_order;

        // Oblique free-wave perturbation, snapped to exact integer Fourier
        // modes so the seed is periodic on the box.
        const Real kx_p_want = crescent_kx_ratio * kc;
        const Real ky_p_want = crescent_ky_ratio * kc;
        const int pmode_x = std::max(1, int(std::round(kx_p_want / dk_x)));
        const int pmode_y = std::max(1, int(std::round(ky_p_want / dk_y)));
        const Real kx_p = pmode_x * dk_x;
        const Real ky_p = pmode_y * dk_y;
        const Real Kp = std::sqrt(kx_p * kx_p + ky_p * ky_p);
        const Real eta_pert_amp = crescent_eps * a;
        const Real phi_pert_amp = crescent_eps * a * std::sqrt(g / Kp);

        amrex::Print() << "  carrier: mode_x=" << wave_mode_x << ", k=" << kc
                       << ", ak=" << ak << ", a=" << a
                       << ", Stokes order=" << order << std::endl;
        amrex::Print() << "  perturbation: modes (" << pmode_x << ", " << pmode_y
                       << "), (kx,ky)=(" << kx_p << ", " << ky_p
                       << "), |k'|=" << Kp << ", eps=" << crescent_eps << std::endl;
        amrex::Print() << "  eta_pert_amp=" << eta_pert_amp
                       << ", phi_pert_amp=" << phi_pert_amp << std::endl;
        if (std::abs(kx_p - kx_p_want) > 1.0e-6 ||
            std::abs(ky_p - ky_p_want) > 1.0e-6) {
            amrex::Print() << "  WARNING: perturbation wavenumbers snapped to grid modes"
                           << " (requested kx=" << kx_p_want << ", ky=" << ky_p_want
                           << "); adjust the domain so kx/ky are integer multiples of"
                           << " dk_x/dk_y to reproduce the reference exactly." << std::endl;
        }

        for (MFIter mfi(eta_o, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto eta_o_arr = eta_o.array(mfi);
            auto phi_o_arr = phi_o.array(mfi);

            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx{
                geom.CellSize()[0], geom.CellSize()[1], geom.CellSize()[2]};
            const auto& domain = geom.Domain();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const Real x = (i - domain.smallEnd(0)) * dx[0];
                const Real y = (j - domain.smallEnd(1)) * dx[1];
                const Real theta = kc * x;

                const Real ak2 = ak * ak;
                const Real ak3 = ak2 * ak;
                const Real ak4 = ak3 * ak;
                const Real ak5 = ak4 * ak;

                // Carrier surface elevation (5th-order Stokes; terms gated by order).
                Real etac = ak * std::cos(theta);
                if (order >= 2) etac += 0.5 * ak2 * std::cos(2.0 * theta);
                if (order >= 3) etac += (3.0 / 8.0) * ak3 *
                                        (std::cos(3.0 * theta) - std::cos(theta));
                if (order >= 4) etac += (1.0 / 3.0) * ak4 *
                                        (std::cos(2.0 * theta) + std::cos(4.0 * theta));
                if (order >= 5) etac += (1.0 / 384.0) * ak5 *
                                        (-422.0 * std::cos(theta) +
                                          297.0 * std::cos(3.0 * theta) +
                                          125.0 * std::cos(5.0 * theta));
                etac *= (1.0 / kc);

                // Carrier surface potential (depth factors use the carrier-only etac).
                const Real ekh = std::exp(kc * etac);
                Real phic = ak * ekh * std::sin(theta);
                if (order >= 3) phic += -0.5 * ak3 * ekh * std::sin(theta);
                if (order >= 4) phic += 0.5 * ak4 * std::exp(2.0 * kc * etac) *
                                        std::sin(2.0 * theta);
                if (order >= 5) phic += ak5 *
                                        (-37.0 / 24.0 * ekh * std::sin(theta) +
                                          1.0 / 12.0 * std::exp(3.0 * kc * etac) *
                                          std::sin(3.0 * theta));
                phic *= sqrt_g_over_k3;

                // Oblique perturbation (potential added first; elevation subtracted last).
                phi_o_arr(i, j, k) = phic + phi_pert_amp *
                                     std::cos(kx_p * x) * std::cos(ky_p * y) *
                                     std::exp(Kp * etac);
                eta_o_arr(i, j, k) = etac - eta_pert_amp *
                                     std::sin(kx_p * x) * std::cos(ky_p * y);
            });
        }
    }

    ParallelCopy(eta, eta_o, 0, 0, 1);
    ParallelCopy(phi, phi_o, 0, 0, 1);

    amrex::Print() << "=== END CRESCENT ===" << std::endl << std::endl;
}

void initializeCoreWaveField(
    crescent::CoreSolverState& core,
    crescent::shared::LevelMeshView level_mesh,
    int finest_level)
{
    const auto& initialization = core.wave.initialization;
    amrex::Print() << "\n=== Wave Initialization: "
                   << initialization.wave_type << " ===" << std::endl;

    if (initialization.wave_type == "stokes") {
        initializeCoreStokesWave(core, level_mesh);
    } else if (initialization.wave_type == "crescent") {
        initializeCoreCrescentWave(core, level_mesh);
    } else if (initialization.wave_type == "jonswap") {
        initializeCoreJonswapWaves(core, level_mesh);
        configureCoreJonswapTimeStep(core.runtime, core.wave, finest_level);
    } else {
        amrex::Abort("Unknown wave_type: " + initialization.wave_type +
                     ". Valid options: 'stokes', 'crescent', 'jonswap'");
    }

    core.solution.saveSurfaceState(0);
}

} // namespace crescent::init

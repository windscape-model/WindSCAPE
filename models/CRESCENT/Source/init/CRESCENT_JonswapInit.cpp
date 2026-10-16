#include "CRESCENT_Init.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_JonswapSpectrum.H"
#include "CRESCENT_Spectral.H"

#include <AMReX_GpuComplex.H>
#include <AMReX_Loop.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>

using namespace amrex;
using namespace crescent::constants;

namespace {

struct JonswapPeakData
{
    Real alpha;
    Real peak_frequency;
    Real peak_wavelength;
};

JonswapPeakData computeJonswapPeakData(const Real wind_speed_10m, const Real fetch)
{
    JonswapPeakData peak_data;
    peak_data.alpha = 0.076 * std::pow((wind_speed_10m*wind_speed_10m) / fetch / GRAVITY, 0.22);
    peak_data.peak_frequency = 22.0 * std::pow((GRAVITY*GRAVITY / wind_speed_10m / fetch),(1.0/3.0));
    peak_data.peak_wavelength = PI2 * GRAVITY / (peak_data.peak_frequency*peak_data.peak_frequency);
    return peak_data;
}

void initializeJonswapSpectrum(
    cMultiFab& eta1,
    cMultiFab& phi1,
    const int nx,
    const int ny,
    const Real dk_x,
    const Real dk_y,
    const Real inverse_gravity,
    const Real alpha,
    const Real peak_frequency,
    const Real initial_angle,
    const Real angle_spread,
    const int directional_exponent,
    const int use_directional_spectrum,
    const Real amplitude_scale)
{
    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    const Real directional_normalization = use_directional_spectrum == 1
        ? crescent::init::jonswap::cosinePowerNormalization(
              angle_spread, directional_exponent)
        : 0.0;
    const Real wrapped_initial_angle = std::remainder(initial_angle, PI2);
    const crescent::init::jonswap::SpectrumParameters parameters{
        nx,
        dk_x,
        dk_y,
        inverse_gravity,
        alpha,
        peak_frequency,
        wrapped_initial_angle,
        angle_spread,
        directional_normalization,
        amplitude_scale,
        directional_exponent,
        use_directional_spectrum};
    const int positive_ky_count = (ny + 1) / 2;

    // The phases are drawn on the host, in one fixed order, so that CPU and
    // GPU builds make the same sea from the same seed (the device generator
    // would make a different sea). The order is that of the former CPU loop:
    // the half-plane is drawn once per tile of the CPU tiling of the single
    // root box (tile size 1024000 x 8 x 8, so max(ny/8, 1) passes), and the
    // last pass remains. That repetition defines the seas of the existing
    // decks, so it stays.
    const IntVect cpu_tile_size(AMREX_D_DECL(1024000, 8, 8));
    for (MFIter mfi(eta1); mfi.isValid(); ++mfi)
    {
        const Box& root = mfi.validbox();
        int passes = 1;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            passes *= std::max(root.length(d) / cpu_tile_size[d], 1);
        }
        Box bx = root;
        bx.setRange(0, 0, nx/2 + 1);
        bx.setRange(1, 0, positive_ky_count);

        eta1[mfi].setVal<RunOn::Host>(GpuComplex<Real>(0.0, 0.0));
        phi1[mfi].setVal<RunOn::Host>(GpuComplex<Real>(0.0, 0.0));
        auto eta1_arr = eta1.array(mfi);
        auto phi1_arr = phi1.array(mfi);

        for (int pass = 0; pass < passes; ++pass) {
            LoopOnCpu(bx, [&] (int i, int j, int z) noexcept
            {
                const Real positive_phase = amrex::Random() * PI2;
                const Real negative_phase = amrex::Random() * PI2;
                Real positive_opposite_phase = positive_phase;
                Real negative_opposite_phase = negative_phase;
                const bool is_zero_mode = i == 0 && j == 0;
                const bool is_x_nyquist =
                    parameters.nx % 2 == 0 && i == parameters.nx / 2;
                if (!is_zero_mode && !is_x_nyquist) {
                    const Real k_x = i * parameters.dk_x;
                    const Real k_y = j * parameters.dk_y;
                    const bool positive_forward =
                        crescent::init::jonswap::hasDirectionalSupport(
                            k_x, k_y, parameters);
                    const bool positive_opposite =
                        crescent::init::jonswap::hasDirectionalSupport(
                            -k_x, -k_y, parameters);
                    if (positive_forward && positive_opposite) {
                        positive_opposite_phase = amrex::Random() * PI2;
                    }

                    if (i > 0 && j > 0) {
                        const bool negative_forward =
                            crescent::init::jonswap::hasDirectionalSupport(
                                k_x, -k_y, parameters);
                        const bool negative_opposite =
                            crescent::init::jonswap::hasDirectionalSupport(
                                -k_x, k_y, parameters);
                        if (negative_forward && negative_opposite) {
                            negative_opposite_phase = amrex::Random() * PI2;
                        }
                    }
                }
                const auto modes = crescent::init::jonswap::makeModePair(
                    i,
                    j,
                    parameters,
                    positive_phase,
                    negative_phase,
                    positive_opposite_phase,
                    negative_opposite_phase);

                eta1_arr(i, j, z) = modes.positive_ky.eta;
                phi1_arr(i, j, z) = modes.positive_ky.phi;
                if (j > 0) {
                    eta1_arr(i, ny - j, z) = modes.negative_ky.eta;
                    phi1_arr(i, ny - j, z) = modes.negative_ky.phi;
                }
            });
        }
    }
}

} // namespace

namespace crescent::init {

void generateJonswapWaves(MultiFab& eta_saved,
                          MultiFab& phi_saved,
                          std::unique_ptr<FFT::R2C<Real>>& fft,
                          const GeometryData& geom,
                          const BoxArray& cba_root,
                          const DistributionMapping& cdmap_root,
                          const Real dk_x,
                          const Real dk_y,
                          const Real wind_speed_10m,
                          const Real fetch,
                          const Real inverse_gravity,
                          const Real initial_angle,
                          const Real angle_spread,
                          const int directional_exponent,
                          const int use_directional_spectrum,
                          const Real jonswap_amplitude_scale,
                          Real& alpha,
                          Real& peak_frequency,
                          Real& peak_wavelength)
{
    const JonswapPeakData peak_data = computeJonswapPeakData(wind_speed_10m, fetch);
    if (!(peak_data.alpha > 0.0) || !std::isfinite(peak_data.alpha) ||
        !(peak_data.peak_frequency > 0.0) ||
        !std::isfinite(peak_data.peak_frequency) ||
        !(peak_data.peak_wavelength > 0.0) ||
        !std::isfinite(peak_data.peak_wavelength)) {
        amrex::Abort(
            "JONSWAP wind_speed_10m/fetch produced non-finite peak parameters");
    }
    alpha = peak_data.alpha;
    peak_frequency = peak_data.peak_frequency;
    peak_wavelength = peak_data.peak_wavelength;

    auto const &[cba, cdm] = fft->getSpectralDataLayout();
    // Host-visible on every build: initializeJonswapSpectrum fills them on
    // the host.
    cMultiFab eta1(cba_root, cdmap_root, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    cMultiFab phi1(cba_root, cdmap_root, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));

    cMultiFab eta0_c(cba, cdm, 1, 0);
    cMultiFab phi0_c(cba, cdm, 1, 0);

    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);

    initializeJonswapSpectrum(
        eta1,
        phi1,
        nx,
        ny,
        dk_x,
        dk_y,
        inverse_gravity,
        alpha,
        peak_frequency,
        initial_angle,
        angle_spread,
        directional_exponent,
        use_directional_spectrum,
        jonswap_amplitude_scale);

    ParallelDescriptor::Barrier();

    ParallelCopy(eta0_c, eta1, 0, 0, 1);
    ParallelCopy(phi0_c, phi1, 0, 0, 1);

    crescent::spectral::backwardFFT(fft, eta0_c, eta_saved);
    crescent::spectral::backwardFFT(fft, phi0_c, phi_saved);

    // Init-time local scratch; hot-path callers use the persistent
    // SpectralWorkspace dealiasing buffer instead.
    amrex::MultiFab dealiasing_scratch(
        eta_saved.boxArray(),
        eta_saved.DistributionMap(),
        1,
        eta_saved.nGrow());
    crescent::spectral::dealiasing(eta_saved, fft, geom, dealiasing_scratch);
    crescent::spectral::dealiasing(phi_saved, fft, geom, dealiasing_scratch);
}

} // namespace crescent::init

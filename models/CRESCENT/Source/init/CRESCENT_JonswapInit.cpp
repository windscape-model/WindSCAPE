#include "CRESCENT_Init.H"
#include "CRESCENT_Constants.H"
#include "CRESCENT_Spectral.H"

#include <AMReX_GpuComplex.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

using namespace amrex;
using namespace crescent::constants;

namespace {

struct JonswapPeakData
{
    Real alpha;
    Real peak_frequency;
    Real peak_wavelength;
};

struct SpectralKyEnergy
{
    Real zero = 0.0;
    Real positive = 0.0;
    Real negative = 0.0;
    Real nyquist = 0.0;
};

JonswapPeakData computeJonswapPeakData(const Real wind_speed_10m, const Real fetch)
{
    JonswapPeakData peak_data;
    peak_data.alpha = 0.076 * std::pow((wind_speed_10m*wind_speed_10m) / fetch / GRAVITY, 0.22);
    peak_data.peak_frequency = 22.0 * std::pow((GRAVITY*GRAVITY / wind_speed_10m / fetch),(1.0/3.0));
    peak_data.peak_wavelength = PI2 * GRAVITY / (peak_data.peak_frequency*peak_data.peak_frequency);
    return peak_data;
}

void accumulateSpectralKyEnergy(const cMultiFab& spectral_field,
                                const int ny,
                                const bool include_nyquist,
                                SpectralKyEnergy& energy)
{
    for (MFIter mfi(spectral_field); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto arr = spectral_field.const_array(mfi);

        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            {
                int j = 0;
                GpuComplex<Real> val = arr(i, j, 0);
                energy.zero += val.real() * val.real() + val.imag() * val.imag();
            }

            for (int j = 1; j < ny/2; ++j) {
                GpuComplex<Real> val = arr(i, j, 0);
                energy.positive += val.real() * val.real() + val.imag() * val.imag();
            }

            if (include_nyquist) {
                int j = ny/2;
                GpuComplex<Real> val = arr(i, j, 0);
                energy.nyquist += val.real() * val.real() + val.imag() * val.imag();
            }

            for (int j = ny/2 + 1; j < ny; ++j) {
                GpuComplex<Real> val = arr(i, j, 0);
                energy.negative += val.real() * val.real() + val.imag() * val.imag();
            }
        }
    }
}

void checkSpectralSymmetry(
    const cMultiFab& spectral_field,
    int ny,
    const std::string& field_name,
    int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    SpectralKyEnergy energy;
    accumulateSpectralKyEnergy(spectral_field, ny, true, energy);

    ParallelDescriptor::ReduceRealSum(energy.zero);
    ParallelDescriptor::ReduceRealSum(energy.positive);
    ParallelDescriptor::ReduceRealSum(energy.negative);
    ParallelDescriptor::ReduceRealSum(energy.nyquist);

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    Print() << "\n=== Spectral Energy Symmetry Check: " << field_name << " ===" << std::endl;
    Print() << "  Energy at k_y = 0:        " << energy.zero << std::endl;
    Print() << "  Energy in +k_y region:    " << energy.positive << std::endl;
    Print() << "  Energy in -k_y region:    " << energy.negative << std::endl;
    Print() << "  Energy at k_y = Nyquist:  " << energy.nyquist << std::endl;
    Print() << "  Ratio (pos/neg):          "
            << (energy.negative > 0 ? energy.positive / energy.negative : 0.0) << std::endl;

    Real total_energy = energy.zero + energy.positive + energy.negative + energy.nyquist;
    Print() << "  Total spectral energy:    " << total_energy << std::endl;

    if (std::abs(energy.positive - energy.negative) / total_energy > 0.01) {
        Print() << "  WARNING: Energy asymmetry > 1% detected!" << std::endl;
        Print() << "  This indicates directional spreading is not symmetric." << std::endl;
    } else {
        Print() << "  OK: Energy is symmetric within 1%" << std::endl;
    }
    Print() << "============================================\n" << std::endl;
}

void initializeJonswapSpectrum(cMultiFab& eta1,
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
                               int verbosity)
{
    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    for (MFIter mfi(eta1, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        bx.setRange(0, 0, nx/2 + 1);
        bx.setRange(1, 0, ny/2);

        auto eta1_arr = eta1.array(mfi);
        auto phi1_arr = phi1.array(mfi);

        const Real dk_x_ = dk_x;
        const Real dk_y_ = dk_y;

        ParallelForRNG(bx, [=] AMREX_GPU_DEVICE(int i, int j, int z, const RandomEngine& engine) noexcept
        {
            if (i > 0 || j > 0) {
                Real k_x = i * dk_x_;
                Real k_y = j * dk_y_;
                Real k_mag = std::sqrt(k_x*k_x + k_y*k_y);

                Real theta_pos = amrex::Random(engine) * PI2;
                Real theta_neg = amrex::Random(engine) * PI2;

                Real omega = sqrt(k_mag / inverse_gravity);
                Real sigma;
                if (omega < peak_frequency) {
                    sigma = 0.07;
                } else {
                    sigma = 0.09;
                }

                Real peak_enhancement = std::exp(-std::pow((omega - peak_frequency),2)/ 2.0 / std::pow(sigma,2) / std::pow(peak_frequency,2));
                Real jonswap_base = alpha / (inverse_gravity*inverse_gravity) / std::pow(omega,5) * std::exp(-1.25*std::pow((omega/peak_frequency),-4.0)) * std::pow(3.3, peak_enhancement);
                jonswap_base = jonswap_base / (inverse_gravity*inverse_gravity) / 2.0 / std::pow(omega, 3);

                Real directional_factor_pos, directional_factor_neg;

                if (angle_spread < 1e-12) {
                    directional_factor_pos = (j == 0) ? 1.0 : 0.0;
                    directional_factor_neg = (j == 0) ? 1.0 : 0.0;
                } else {
                    Real wave_angle = (j > 0) ? std::atan2(k_y, k_x) : 0.0;
                    if (wave_angle < 0.0) wave_angle += PI2;
                    Real angle_diff = wave_angle - initial_angle;
                    while (angle_diff > PI) angle_diff -= PI2;
                    while (angle_diff < -PI) angle_diff += PI2;

                    Real abs_angle_diff = std::abs(angle_diff);

                    if (abs_angle_diff <= 0.5*angle_spread) {
                        Real directional_factor = 2.0/angle_spread * std::pow(std::cos(PI2/2.0*abs_angle_diff/angle_spread), directional_exponent);
                        directional_factor_pos = directional_factor;
                        directional_factor_neg = directional_factor;
                    } else {
                        directional_factor_pos = 0.0;
                        directional_factor_neg = 0.0;
                    }
                }

                Real spectrum_value_pos = jonswap_base * directional_factor_pos;
                Real spectrum_value_neg = jonswap_base * directional_factor_neg;

                Real wave_amplitude_pos = sqrt(2.0 * spectrum_value_pos * dk_x_ * dk_y_);
                Real wave_amplitude_neg = sqrt(2.0 * spectrum_value_neg * dk_x_ * dk_y_);

                wave_amplitude_pos *= 0.25;
                wave_amplitude_neg *= 0.25;
                if (i == 0 || j == 0) {
                    wave_amplitude_pos *= sqrt(2.0);
                    wave_amplitude_neg *= sqrt(2.0);
                }

                Real phi_scale_pos = wave_amplitude_pos / inverse_gravity / omega;
                Real phi_scale_neg = wave_amplitude_neg / inverse_gravity / omega;

                Real cos_theta_pos = std::cos(theta_pos);
                Real sin_theta_pos = std::sin(theta_pos);
                Real cos_theta_neg = std::cos(theta_neg);
                Real sin_theta_neg = std::sin(theta_neg);

                GpuComplex<Real> eta_pos = GpuComplex<Real>(
                    wave_amplitude_pos * cos_theta_pos,
                    wave_amplitude_pos * sin_theta_pos);
                GpuComplex<Real> phi_pos = GpuComplex<Real>(
                    +phi_scale_pos * sin_theta_pos,
                    -phi_scale_pos * cos_theta_pos);

                GpuComplex<Real> eta_neg = GpuComplex<Real>(
                    wave_amplitude_neg * cos_theta_neg,
                    wave_amplitude_neg * sin_theta_neg);
                GpuComplex<Real> phi_neg = GpuComplex<Real>(
                    +phi_scale_neg * sin_theta_neg,
                    -phi_scale_neg * cos_theta_neg);

                eta1_arr(i, j, z) = eta_pos;
                phi1_arr(i, j, z) = phi_pos;

                if (j > 0) {
                    eta1_arr(i, ny - j, z) = eta_neg;
                    phi1_arr(i, ny - j, z) = phi_neg;
                }

                if (verbosity >= 2 && i == 10 && j == 1) {
                    printf("Mode (%d,%d): theta_pos=%.6f, theta_neg=%.6f rad, amp_pos=%.6e, amp_neg=%.6e, eta_pos=(%.6e,%.6e), eta_neg=(%.6e,%.6e)\n",
                           i, j, theta_pos, theta_neg,
                           wave_amplitude_pos, wave_amplitude_neg,
                           eta_pos.real(), eta_pos.imag(),
                           eta_neg.real(), eta_neg.imag());
                }
            }
        });
    }
}

void printDirectionalEnergyDistribution(
    const cMultiFab& eta1,
    const int ny,
    int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    Print() << "\n=== JONSWAP Directional Energy Distribution ===" << std::endl;

    SpectralKyEnergy energy;
    accumulateSpectralKyEnergy(eta1, ny, false, energy);

    ParallelDescriptor::ReduceRealSum(energy.zero);
    ParallelDescriptor::ReduceRealSum(energy.positive);
    ParallelDescriptor::ReduceRealSum(energy.negative);

    Real total_energy = energy.zero + energy.positive + energy.negative;

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    Print() << "  Energy at ky = 0:          " << energy.zero
            << " (" << 100.0 * energy.zero / total_energy << "%)" << std::endl;
    Print() << "  Energy in +ky region:      " << energy.positive
            << " (" << 100.0 * energy.positive / total_energy << "%)" << std::endl;
    Print() << "  Energy in -ky region:      " << energy.negative
            << " (" << 100.0 * energy.negative / total_energy << "%)" << std::endl;
    Print() << "  Ratio (+ky / -ky):         "
            << (energy.negative > 0 ? energy.positive / energy.negative : 0.0) << std::endl;
    Print() << "  Total spectral energy:     " << total_energy << std::endl;

    if (std::abs(energy.positive - energy.negative) / total_energy < 0.01) {
        Print() << "  ✓ Energy is symmetric between ±ky within 1%" << std::endl;
    } else {
        Print() << "  ⚠ Energy asymmetry > 1% (expected for initial_angle ≠ 0)" << std::endl;
    }
    Print() << "============================================\n" << std::endl;
}

void enforceYSymmetryIfEnabled(MultiFab& eta_saved, MultiFab& phi_saved, const GeometryData& geom)
{
    const bool enforce_y_symmetry = false;
    if (!enforce_y_symmetry) {
        return;
    }

    const int ny = geom.Domain().length(1);

    amrex::Print() << "\n=== Enforcing y-symmetry in physical space ===" << std::endl;
    amrex::Print() << "  Domain: ny=" << ny << std::endl;

    MultiFab eta_temp(eta_saved.boxArray(), eta_saved.DistributionMap(), 1, 0);
    MultiFab phi_temp(phi_saved.boxArray(), phi_saved.DistributionMap(), 1, 0);

    for (MFIter mfi(eta_saved, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto eta_arr = eta_saved.const_array(mfi);
        auto phi_arr = phi_saved.const_array(mfi);
        auto eta_temp_arr = eta_temp.array(mfi);
        auto phi_temp_arr = phi_temp.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int j_mirror = (ny - j) % ny;

            Real eta_sym = 0.5 * (eta_arr(i, j, k) + eta_arr(i, j_mirror, k));
            Real phi_sym = 0.5 * (phi_arr(i, j, k) + phi_arr(i, j_mirror, k));

            eta_temp_arr(i, j, k) = eta_sym;
            phi_temp_arr(i, j, k) = phi_sym;
        });
    }

    MultiFab::Copy(eta_saved, eta_temp, 0, 0, 1, 0);
    MultiFab::Copy(phi_saved, phi_temp, 0, 0, 1, 0);

    amrex::Print() << "  After symmetrization:" << std::endl;
    amrex::Print() << "    eta min/max: " << eta_saved.min(0) << " / " << eta_saved.max(0) << std::endl;
    amrex::Print() << "    phi min/max: " << phi_saved.min(0) << " / " << phi_saved.max(0) << std::endl;

    Real max_eta_asym = 0.0;
    Real max_phi_asym = 0.0;
    for (MFIter mfi(eta_saved); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto eta_arr = eta_saved.const_array(mfi);
        auto phi_arr = phi_saved.const_array(mfi);

        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                int j_mirror = (ny - j) % ny;
                Real eta_diff = std::abs(eta_arr(i, j, 0) - eta_arr(i, j_mirror, 0));
                Real phi_diff = std::abs(phi_arr(i, j, 0) - phi_arr(i, j_mirror, 0));
                max_eta_asym = std::max(max_eta_asym, eta_diff);
                max_phi_asym = std::max(max_phi_asym, phi_diff);
            }
        }
    }
    ParallelDescriptor::ReduceRealMax(max_eta_asym);
    ParallelDescriptor::ReduceRealMax(max_phi_asym);

    amrex::Print() << "  Symmetry verification:" << std::endl;
    amrex::Print() << "    Max |eta(x,y) - eta(x,-y)|: " << max_eta_asym << std::endl;
    amrex::Print() << "    Max |phi(x,y) - phi(x,-y)|: " << max_phi_asym << std::endl;
    if (max_eta_asym < 1e-12 && max_phi_asym < 1e-12) {
        amrex::Print() << "    ✓ Perfect symmetry (< 1e-12)" << std::endl;
    } else if (max_eta_asym < 1e-10) {
        amrex::Print() << "    ✓ Excellent symmetry (< 1e-10)" << std::endl;
    } else {
        amrex::Print() << "    ⚠ Symmetry error > 1e-10" << std::endl;
    }
    amrex::Print() << "===============================================\n" << std::endl;
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
                          Real& alpha,
                          Real& peak_frequency,
                          Real& peak_wavelength,
                          int verbosity)
{
    const JonswapPeakData peak_data = computeJonswapPeakData(wind_speed_10m, fetch);
    alpha = peak_data.alpha;
    peak_frequency = peak_data.peak_frequency;
    peak_wavelength = peak_data.peak_wavelength;

    auto const &[cba, cdm] = fft->getSpectralDataLayout();
    cMultiFab eta1(cba_root, cdmap_root, 1, 0);
    cMultiFab phi1(cba_root, cdmap_root, 1, 0);

    eta1.setVal(0.0);
    phi1.setVal(0.0);

    cMultiFab eta0_c(cba, cdm, 1, 0);
    cMultiFab phi0_c(cba, cdm, 1, 0);

    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);

    initializeJonswapSpectrum(
        eta1, phi1, nx, ny, dk_x, dk_y, inverse_gravity,
        alpha,
        peak_frequency,
        initial_angle,
        angle_spread,
        directional_exponent,
        verbosity);

    ParallelDescriptor::Barrier();

    checkSpectralSymmetry(eta1, ny, "eta (before inverse FFT)", verbosity);
    checkSpectralSymmetry(phi1, ny, "phi (before inverse FFT)", verbosity);
    printDirectionalEnergyDistribution(eta1, ny, verbosity);

    ParallelCopy(eta0_c, eta1, 0, 0, 1);
    ParallelCopy(phi0_c, phi1, 0, 0, 1);

    crescent::spectral::backwardFFT(fft, geom, eta0_c, eta_saved);
    crescent::spectral::backwardFFT(fft, geom, phi0_c, phi_saved);

    // Init-time local scratch; hot-path callers use the persistent
    // SpectralWorkspace dealiasing buffer instead.
    amrex::MultiFab dealiasing_scratch(
        eta_saved.boxArray(),
        eta_saved.DistributionMap(),
        1,
        eta_saved.nGrow());
    crescent::spectral::dealiasing(eta_saved, fft, geom, dealiasing_scratch);
    crescent::spectral::dealiasing(phi_saved, fft, geom, dealiasing_scratch);

    enforceYSymmetryIfEnabled(eta_saved, phi_saved, geom);
}

} // namespace crescent::init

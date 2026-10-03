#pragma once

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace crescent::spectral::filter_kernels {

struct LowPassFilterParameters
{
    int ny;
    amrex::Real dk_x;
    amrex::Real dk_y;
    amrex::Real normalization_factor;
    amrex::Real k_peak;
    amrex::Real beta1;
    amrex::Real beta2;
    bool zero_mean;
    bool require_positive_peak;
};

struct CombinedSpectralFilterParameters
{
    int ny;
    amrex::Real dk_x;
    amrex::Real dk_y;
    amrex::Real normalization_factor;
    amrex::Real k_nyquist;
    amrex::Real k_dealiasing;
    amrex::Real k_peak;
    amrex::Real beta1;
    amrex::Real beta2;
    amrex::Real houli_alpha;
    int houli_power;
    bool zero_mean;
    bool require_positive_peak;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int spectralKxIndex(int i, int nx) noexcept
{
    return (i <= nx / 2) ? i : (i - nx);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int spectralKyIndex(int j, int ny) noexcept
{
    return (j <= ny / 2) ? j : (j - ny);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool spectralModeExceedsCutoff(int kx_index,
                               int ky_index,
                               int cutoff_x,
                               int cutoff_y) noexcept
{
    const int abs_kx = (kx_index < 0) ? -kx_index : kx_index;
    const int abs_ky = (ky_index < 0) ? -ky_index : ky_index;
    return abs_kx > cutoff_x || abs_ky > cutoff_y;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real spectralMagnitude(int i,
                              int j,
                              int ny,
                              amrex::Real dk_x,
                              amrex::Real dk_y) noexcept
{
    const amrex::Real k_x = i * dk_x;
    const amrex::Real k_y = spectralKyIndex(j, ny) * dk_y;
    return std::sqrt(k_x * k_x + k_y * k_y);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real spectralMagnitudeFromIndices(int kx_index,
                                         int ky_index,
                                         amrex::Real dk_x,
                                         amrex::Real dk_y) noexcept
{
    const amrex::Real k_x = kx_index * dk_x;
    const amrex::Real k_y = ky_index * dk_y;
    return std::sqrt(k_x * k_x + k_y * k_y);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real lowPassWeight(amrex::Real k_mag,
                          amrex::Real k_peak,
                          amrex::Real beta1,
                          amrex::Real beta2,
                          bool require_positive_peak) noexcept
{
    const bool apply_low_pass = require_positive_peak
        ? (k_peak > 0.0 && k_mag > k_peak)
        : (k_mag > k_peak);

    if (!apply_low_pass) {
        return 1.0;
    }

    const amrex::Real exponent =
        std::pow(std::abs(k_mag / beta1 / k_peak), beta2);
    return std::exp(-exponent);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real lowPassWeight(amrex::Real k_mag,
                          const LowPassFilterParameters& params) noexcept
{
    return lowPassWeight(
        k_mag,
        params.k_peak,
        params.beta1,
        params.beta2,
        params.require_positive_peak);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real houLiWeight(amrex::Real k_mag,
                        amrex::Real k_nyquist,
                        amrex::Real houli_alpha,
                        int houli_power) noexcept
{
    return std::exp(-houli_alpha * std::pow(k_mag / k_nyquist, houli_power));
}

template <typename Complex>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void applyLowPassZeroMean(int i,
                          int j,
                          Complex& fft_value,
                          const LowPassFilterParameters& params) noexcept
{
    const amrex::Real k_mag =
        spectralMagnitude(i, j, params.ny, params.dk_x, params.dk_y);

    const bool apply_low_pass = params.require_positive_peak
        ? (params.k_peak > 0.0 && k_mag > params.k_peak)
        : (k_mag > params.k_peak);

    if (apply_low_pass) {
        const amrex::Real sigma = lowPassWeight(k_mag, params);
        fft_value = fft_value * sigma / params.normalization_factor;
    } else if (params.zero_mean && i == 0 && j == 0) {
        fft_value = Complex(0.0, 0.0);
    } else {
        fft_value = fft_value / params.normalization_factor;
    }
}

template <typename Complex>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void applyCombinedHouLiDealiasLowPassZeroMean(
    int i,
    int j,
    Complex& fft_value,
    const CombinedSpectralFilterParameters& params) noexcept
{
    if (params.zero_mean && i == 0 && j == 0) {
        fft_value = Complex(0.0, 0.0);
        return;
    }

    const amrex::Real k_mag =
        spectralMagnitude(i, j, params.ny, params.dk_x, params.dk_y);
    amrex::Real sigma = houLiWeight(
        k_mag,
        params.k_nyquist,
        params.houli_alpha,
        params.houli_power);

    if (k_mag > params.k_dealiasing) {
        sigma = 0.0;
    }

    const LowPassFilterParameters low_pass_params{
        params.ny,
        params.dk_x,
        params.dk_y,
        params.normalization_factor,
        params.k_peak,
        params.beta1,
        params.beta2,
        false,
        params.require_positive_peak
    };
    sigma *= lowPassWeight(k_mag, low_pass_params);

    if (params.zero_mean && i == 0 && j == 0) {
        fft_value = Complex(0.0, 0.0);
    } else {
        fft_value = fft_value * sigma / params.normalization_factor;
    }
}

} // namespace crescent::spectral::filter_kernels

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
    int nx;
    int ny;
    amrex::Real dk_x;
    amrex::Real dk_y;
    amrex::Real normalization_factor;
    amrex::Real k_nyquist;
    // Two-thirds dealiasing as a per-axis (square) mode-index cutoff, the
    // same rule as the parent's dealiasing and DNO recursion.
    int cutoff_x;
    int cutoff_y;
    amrex::Real k_peak;
    amrex::Real beta1;
    amrex::Real beta2;
    amrex::Real houli_alpha;
    int houli_power;
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

// Parent low-pass weight of mode (i, j): the low-pass weight above the
// peak, 1 below it.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real lowPassModeWeight(int i,
                              int j,
                              const LowPassFilterParameters& params) noexcept
{
    const amrex::Real k_mag =
        spectralMagnitude(i, j, params.ny, params.dk_x, params.dk_y);

    const bool apply_low_pass = params.require_positive_peak
        ? (params.k_peak > 0.0 && k_mag > params.k_peak)
        : (k_mag > params.k_peak);

    return apply_low_pass ? lowPassWeight(k_mag, params) : amrex::Real(1.0);
}

// Child filter weight of mode (i, j): Hou-Li roll-off times the square
// two-thirds cut times the low-pass weight. The k = 0 mode is kept (a child
// window has a nonzero mean).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real combinedFilterWeight(
    int i,
    int j,
    const CombinedSpectralFilterParameters& params) noexcept
{
    const amrex::Real k_mag =
        spectralMagnitude(i, j, params.ny, params.dk_x, params.dk_y);
    amrex::Real sigma = houLiWeight(
        k_mag,
        params.k_nyquist,
        params.houli_alpha,
        params.houli_power);

    if (spectralModeExceedsCutoff(
            spectralKxIndex(i, params.nx),
            spectralKyIndex(j, params.ny),
            params.cutoff_x,
            params.cutoff_y)) {
        sigma = 0.0;
    }

    sigma *= lowPassWeight(
        k_mag,
        params.k_peak,
        params.beta1,
        params.beta2,
        params.require_positive_peak);
    return sigma;
}

} // namespace crescent::spectral::filter_kernels

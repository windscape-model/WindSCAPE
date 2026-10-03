#ifndef CRESCENT_RUNTIME_DIAGNOSTIC_METRICS_H
#define CRESCENT_RUNTIME_DIAGNOSTIC_METRICS_H

#include <AMReX_Math.H>
#include <AMReX_Real.H>

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace crescent::solver::runtime_diagnostics::detail {

inline std::string runtimeOutputDirectory(const std::string& plot_file)
{
    const std::size_t slash = plot_file.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string(".") : plot_file.substr(0, slash);
}

struct LineSpectrumMetrics
{
    int n = 0;
    amrex::Real mode0_energy = 0.0;
    amrex::Real nonzero_energy = 0.0;
    amrex::Real total_energy = 0.0;
    amrex::Real nonzero_fraction = 0.0;
    amrex::Real mean_abs_mode = 0.0;
    amrex::Real peak_abs_mode_energy = 0.0;
    int peak_abs_mode = 0;
};

inline LineSpectrumMetrics summarizeLineSpectrum(const std::vector<amrex::Real>& samples)
{
    LineSpectrumMetrics metrics;
    metrics.n = static_cast<int>(samples.size());
    if (metrics.n <= 1) {
        return metrics;
    }

    std::vector<amrex::Real> abs_mode_energy(metrics.n / 2 + 1, 0.0);
    const amrex::Real norm_factor = amrex::Real(1.0) /
                                    (static_cast<amrex::Real>(metrics.n) *
                                     static_cast<amrex::Real>(metrics.n));

    for (int mode = 0; mode < metrics.n; ++mode) {
        std::complex<amrex::Real> coeff(0.0, 0.0);
        for (int idx = 0; idx < metrics.n; ++idx) {
            const amrex::Real angle =
                -amrex::Math::pi<amrex::Real>() * amrex::Real(2.0) *
                static_cast<amrex::Real>(mode * idx) /
                static_cast<amrex::Real>(metrics.n);
            coeff += static_cast<amrex::Real>(samples[idx]) *
                     std::complex<amrex::Real>(std::cos(angle), std::sin(angle));
        }

        const int abs_mode = (mode <= metrics.n / 2) ? mode : (metrics.n - mode);
        abs_mode_energy[abs_mode] += std::norm(coeff) * norm_factor;
    }

    metrics.mode0_energy = abs_mode_energy[0];
    amrex::Real weighted_abs_mode_sum = 0.0;

    for (int abs_mode = 1; abs_mode < static_cast<int>(abs_mode_energy.size()); ++abs_mode) {
        const amrex::Real energy = abs_mode_energy[abs_mode];
        metrics.nonzero_energy += energy;
        weighted_abs_mode_sum += static_cast<amrex::Real>(abs_mode) * energy;
        if (energy > metrics.peak_abs_mode_energy) {
            metrics.peak_abs_mode_energy = energy;
            metrics.peak_abs_mode = abs_mode;
        }
    }

    metrics.total_energy = metrics.mode0_energy + metrics.nonzero_energy;
    metrics.nonzero_fraction =
        (metrics.total_energy > 0.0) ? (metrics.nonzero_energy / metrics.total_energy) : 0.0;
    metrics.mean_abs_mode =
        (metrics.nonzero_energy > 0.0) ? (weighted_abs_mode_sum / metrics.nonzero_energy) : 0.0;
    return metrics;
}

template <typename ArrayType>
std::vector<amrex::Real> extractRowSamples(
    const ArrayType& arr,
    int i_lo,
    int i_hi,
    int j_probe)
{
    std::vector<amrex::Real> row;
    row.reserve(i_hi - i_lo + 1);
    for (int i = i_lo; i <= i_hi; ++i) {
        row.push_back(arr(i, j_probe, 0));
    }
    return row;
}

template <typename ArrayType>
std::vector<amrex::Real> extractColumnSamples(
    const ArrayType& arr,
    int j_lo,
    int j_hi,
    int i_probe)
{
    std::vector<amrex::Real> col;
    col.reserve(j_hi - j_lo + 1);
    for (int j = j_lo; j <= j_hi; ++j) {
        col.push_back(arr(i_probe, j, 0));
    }
    return col;
}

struct JumpMetrics
{
    int n = 0;
    amrex::Real min_jump = 0.0;
    amrex::Real max_jump = 0.0;
    amrex::Real max_abs_jump = 0.0;
    amrex::Real mean_abs_jump = 0.0;
    amrex::Real rms_jump = 0.0;
};

inline JumpMetrics summarizeJumpMetrics(
    const std::vector<amrex::Real>& values0,
    const std::vector<amrex::Real>& values1)
{
    JumpMetrics metrics;
    metrics.n = static_cast<int>(std::min(values0.size(), values1.size()));
    if (metrics.n <= 0) {
        return metrics;
    }

    metrics.min_jump = values1[0] - values0[0];
    metrics.max_jump = metrics.min_jump;
    metrics.max_abs_jump = std::abs(metrics.min_jump);
    amrex::Real sum_abs_jump = 0.0;
    amrex::Real sum_sq_jump = 0.0;

    for (int idx = 0; idx < metrics.n; ++idx) {
        const amrex::Real jump = values1[idx] - values0[idx];
        metrics.min_jump = std::min(metrics.min_jump, jump);
        metrics.max_jump = std::max(metrics.max_jump, jump);
        metrics.max_abs_jump = std::max(metrics.max_abs_jump, std::abs(jump));
        sum_abs_jump += std::abs(jump);
        sum_sq_jump += jump * jump;
    }

    metrics.mean_abs_jump = sum_abs_jump / static_cast<amrex::Real>(metrics.n);
    metrics.rms_jump = std::sqrt(sum_sq_jump / static_cast<amrex::Real>(metrics.n));
    return metrics;
}

struct PatchMetrics
{
    int nx_patch = 0;
    int ny_patch = 0;
    int n = 0;
    amrex::Real min_val = 0.0;
    amrex::Real max_val = 0.0;
    amrex::Real max_abs_val = 0.0;
    amrex::Real rms_val = 0.0;
    amrex::Real max_abs_dx = 0.0;
    amrex::Real max_abs_dy = 0.0;
};

template <typename ArrayType>
bool computePatchMetrics(
    int i_lo_patch,
    int i_hi_patch,
    int j_lo_patch,
    int j_hi_patch,
    const ArrayType& arr,
    PatchMetrics& metrics)
{
    if (i_lo_patch > i_hi_patch || j_lo_patch > j_hi_patch) {
        return false;
    }

    metrics.nx_patch = i_hi_patch - i_lo_patch + 1;
    metrics.ny_patch = j_hi_patch - j_lo_patch + 1;
    metrics.n = metrics.nx_patch * metrics.ny_patch;
    metrics.min_val = arr(i_lo_patch, j_lo_patch, 0);
    metrics.max_val = metrics.min_val;
    metrics.max_abs_val = std::abs(metrics.min_val);
    amrex::Real sum_sq = 0.0;

    for (int j = j_lo_patch; j <= j_hi_patch; ++j) {
        for (int i = i_lo_patch; i <= i_hi_patch; ++i) {
            const amrex::Real val = arr(i, j, 0);
            metrics.min_val = std::min(metrics.min_val, val);
            metrics.max_val = std::max(metrics.max_val, val);
            metrics.max_abs_val = std::max(metrics.max_abs_val, std::abs(val));
            sum_sq += val * val;
            if (i < i_hi_patch) {
                metrics.max_abs_dx = std::max(
                    metrics.max_abs_dx,
                    std::abs(arr(i + 1, j, 0) - val));
            }
            if (j < j_hi_patch) {
                metrics.max_abs_dy = std::max(
                    metrics.max_abs_dy,
                    std::abs(arr(i, j + 1, 0) - val));
            }
        }
    }

    metrics.rms_val = std::sqrt(sum_sq / static_cast<amrex::Real>(metrics.n));
    return true;
}

} // namespace crescent::solver::runtime_diagnostics::detail

#endif // CRESCENT_RUNTIME_DIAGNOSTIC_METRICS_H

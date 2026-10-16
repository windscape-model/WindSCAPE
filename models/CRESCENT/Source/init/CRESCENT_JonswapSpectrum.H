#pragma once

#include <AMReX_GpuComplex.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

#include "CRESCENT_Constants.H"

namespace crescent::init::jonswap {

struct SpectrumParameters
{
    int nx = 0;
    amrex::Real dk_x = 0.0;
    amrex::Real dk_y = 0.0;
    amrex::Real inverse_gravity = 0.0;
    amrex::Real alpha = 0.0;
    amrex::Real peak_frequency = 0.0;
    amrex::Real initial_angle = 0.0;
    amrex::Real angle_spread = 0.0;
    amrex::Real directional_normalization = 0.0;
    amrex::Real amplitude_scale = 1.0;
    int directional_exponent = 0;
    int use_directional_spectrum = 1;
};

struct TravelingMode
{
    amrex::GpuComplex<amrex::Real> eta{0.0, 0.0};
    amrex::GpuComplex<amrex::Real> phi{0.0, 0.0};
};

struct ModePair
{
    TravelingMode positive_ky;
    TravelingMode negative_ky;
};

inline amrex::Real cosinePowerIntegral(const int exponent)
{
    using amrex::Real;
    if (exponent > 128) {
        const Real exponent_real = static_cast<Real>(exponent);
        return std::exp(
            Real(0.5) * std::log(crescent::constants::PI) +
            std::lgamma((exponent_real + Real(1.0)) / Real(2.0)) -
            std::lgamma((exponent_real + Real(2.0)) / Real(2.0)));
    }
    return std::sqrt(crescent::constants::PI) *
           std::tgamma(Real(exponent + 1) / Real(2.0)) /
           std::tgamma(Real(exponent + 2) / Real(2.0));
}

inline amrex::Real cosinePowerNormalization(
    const amrex::Real angle_spread,
    const int exponent)
{
    return crescent::constants::PI /
           (angle_spread * cosinePowerIntegral(exponent));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wrapAngle(amrex::Real angle) noexcept
{
    while (angle > crescent::constants::PI) {
        angle -= crescent::constants::PI2;
    }
    while (angle < -crescent::constants::PI) {
        angle += crescent::constants::PI2;
    }
    return angle;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real directionalSpreading(
    const amrex::Real wave_angle,
    const amrex::Real initial_angle,
    const amrex::Real angle_spread,
    const int exponent,
    const amrex::Real normalization) noexcept
{
    const amrex::Real angle_difference =
        std::abs(wrapAngle(wave_angle - initial_angle));
    if (angle_difference > amrex::Real(0.5) * angle_spread) {
        return 0.0;
    }

    const amrex::Real cosine = std::cos(
        crescent::constants::PI * angle_difference / angle_spread);
    const amrex::Real nonnegative_cosine = cosine > 0.0 ? cosine : 0.0;
    return normalization * std::pow(nonnegative_cosine, exponent);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real frequencySpectrum(
    const amrex::Real omega,
    const amrex::Real inverse_gravity,
    const amrex::Real alpha,
    const amrex::Real peak_frequency) noexcept
{
    const amrex::Real sigma = omega < peak_frequency ? 0.07 : 0.09;
    const amrex::Real normalized_offset =
        (omega - peak_frequency) / (sigma * peak_frequency);
    const amrex::Real peak_enhancement =
        std::exp(-amrex::Real(0.5) * normalized_offset * normalized_offset);
    const amrex::Real gravity_squared =
        amrex::Real(1.0) / (inverse_gravity * inverse_gravity);

    return alpha * gravity_squared / std::pow(omega, 5) *
           std::exp(-amrex::Real(1.25) *
                    std::pow(omega / peak_frequency, -4.0)) *
           std::pow(amrex::Real(3.3), peak_enhancement);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real directionalWavenumberDensity(
    const amrex::Real k_x,
    const amrex::Real k_y,
    const SpectrumParameters& parameters) noexcept
{
    const amrex::Real wavenumber = std::sqrt(k_x * k_x + k_y * k_y);
    const amrex::Real omega = std::sqrt(wavenumber / parameters.inverse_gravity);
    const amrex::Real spreading = directionalSpreading(
        std::atan2(k_y, k_x),
        parameters.initial_angle,
        parameters.angle_spread,
        parameters.directional_exponent,
        parameters.directional_normalization);
    const amrex::Real gravity_squared =
        amrex::Real(1.0) /
        (parameters.inverse_gravity * parameters.inverse_gravity);

    // Polar-to-Cartesian Jacobian:
    // S(kx,ky) = S(omega) D(theta) (d omega / d k) / k.
    return frequencySpectrum(
               omega,
               parameters.inverse_gravity,
               parameters.alpha,
               parameters.peak_frequency) *
           spreading * gravity_squared /
           (amrex::Real(2.0) * std::pow(omega, 3));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool hasDirectionalSupport(
    const amrex::Real k_x,
    const amrex::Real k_y,
    const SpectrumParameters& parameters) noexcept
{
    return parameters.use_directional_spectrum == 1 &&
           directionalSpreading(
               std::atan2(k_y, k_x),
               parameters.initial_angle,
               parameters.angle_spread,
               parameters.directional_exponent,
               parameters.directional_normalization) > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real longCrestedWavenumberDensity(
    const amrex::Real k_x,
    const SpectrumParameters& parameters) noexcept
{
    const amrex::Real omega = std::sqrt(k_x / parameters.inverse_gravity);
    const amrex::Real gravity = amrex::Real(1.0) / parameters.inverse_gravity;

    // S_k = S_omega d omega / d k for deep water.
    return frequencySpectrum(
               omega,
               parameters.inverse_gravity,
               parameters.alpha,
               parameters.peak_frequency) *
           gravity / (amrex::Real(2.0) * omega);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real normalizedCoefficientMagnitude(
    const amrex::Real spectral_density,
    const amrex::Real spectral_cell_measure,
    const amrex::Real amplitude_scale) noexcept
{
    return amrex::Real(0.5) * amplitude_scale *
           std::sqrt(amrex::Real(2.0) * spectral_density * spectral_cell_measure);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
TravelingMode makeTravelingMode(
    const amrex::Real coefficient_magnitude,
    const amrex::Real inverse_gravity,
    const amrex::Real omega,
    const amrex::Real phase) noexcept
{
    const amrex::Real cosine = std::cos(phase);
    const amrex::Real sine = std::sin(phase);
    const amrex::Real potential_scale =
        coefficient_magnitude / (inverse_gravity * omega);

    TravelingMode mode;
    mode.eta = amrex::GpuComplex<amrex::Real>(
        coefficient_magnitude * cosine,
        coefficient_magnitude * sine);
    // A +k traveling wave satisfies phi_k = -i (g / omega) eta_k.
    mode.phi = amrex::GpuComplex<amrex::Real>(
        potential_scale * sine,
        -potential_scale * cosine);
    return mode;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
TravelingMode conjugateMode(const TravelingMode& mode) noexcept
{
    TravelingMode conjugate;
    conjugate.eta = amrex::GpuComplex<amrex::Real>(
        mode.eta.real(), -mode.eta.imag());
    conjugate.phi = amrex::GpuComplex<amrex::Real>(
        mode.phi.real(), -mode.phi.imag());
    return conjugate;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
TravelingMode addModes(
    const TravelingMode& first,
    const TravelingMode& second) noexcept
{
    TravelingMode sum;
    sum.eta = amrex::GpuComplex<amrex::Real>(
        first.eta.real() + second.eta.real(),
        first.eta.imag() + second.eta.imag());
    sum.phi = amrex::GpuComplex<amrex::Real>(
        first.phi.real() + second.phi.real(),
        first.phi.imag() + second.phi.imag());
    return sum;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
TravelingMode makeStoredDirectionalMode(
    const amrex::Real k_x,
    const amrex::Real k_y,
    const amrex::Real omega,
    const SpectrumParameters& parameters,
    const amrex::Real forward_phase,
    const amrex::Real opposite_phase) noexcept
{
    const amrex::Real spectral_cell_measure =
        parameters.dk_x * parameters.dk_y;
    const amrex::Real forward_magnitude = normalizedCoefficientMagnitude(
        directionalWavenumberDensity(k_x, k_y, parameters),
        spectral_cell_measure,
        parameters.amplitude_scale);
    TravelingMode stored = makeTravelingMode(
        forward_magnitude,
        parameters.inverse_gravity,
        omega,
        forward_phase);

    // An R2C coefficient also represents the conjugate wavevector: add the
    // wave travelling toward -k with its own phase.
    const amrex::Real opposite_magnitude = normalizedCoefficientMagnitude(
        directionalWavenumberDensity(-k_x, -k_y, parameters),
        spectral_cell_measure,
        parameters.amplitude_scale);
    const TravelingMode opposite_at_physical_wavevector = makeTravelingMode(
        opposite_magnitude,
        parameters.inverse_gravity,
        omega,
        opposite_phase);
    stored = addModes(
        stored,
        conjugateMode(opposite_at_physical_wavevector));
    return stored;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ModePair makeModePair(
    const int i,
    const int positive_j,
    const SpectrumParameters& parameters,
    const amrex::Real positive_phase,
    const amrex::Real negative_phase,
    const amrex::Real positive_opposite_phase = 0.0,
    const amrex::Real negative_opposite_phase = 0.0) noexcept
{
    ModePair modes;

    const bool is_zero_mode = i == 0 && positive_j == 0;
    const bool is_x_nyquist =
        parameters.nx % 2 == 0 && i == parameters.nx / 2;
    if (is_zero_mode || is_x_nyquist) {
        return modes;
    }

    const amrex::Real k_x = i * parameters.dk_x;
    const amrex::Real k_y = positive_j * parameters.dk_y;
    const amrex::Real wavenumber = std::sqrt(k_x * k_x + k_y * k_y);
    const amrex::Real omega = std::sqrt(wavenumber / parameters.inverse_gravity);

    if (parameters.use_directional_spectrum == 0) {
        if (positive_j != 0) {
            return modes;
        }
        const amrex::Real magnitude = normalizedCoefficientMagnitude(
            longCrestedWavenumberDensity(k_x, parameters),
            parameters.dk_x,
            parameters.amplitude_scale);
        modes.positive_ky = makeTravelingMode(
            magnitude,
            parameters.inverse_gravity,
            omega,
            positive_phase);
        return modes;
    }

    modes.positive_ky = makeStoredDirectionalMode(
        k_x,
        k_y,
        omega,
        parameters,
        positive_phase,
        positive_opposite_phase);

    if (positive_j == 0) {
        return modes;
    }
    if (i == 0) {
        // kx=0 belongs to its own real-FFT conjugacy plane.
        modes.negative_ky = conjugateMode(modes.positive_ky);
        return modes;
    }

    modes.negative_ky = makeStoredDirectionalMode(
        k_x,
        -k_y,
        omega,
        parameters,
        negative_phase,
        negative_opposite_phase);
    return modes;
}

} // namespace crescent::init::jonswap

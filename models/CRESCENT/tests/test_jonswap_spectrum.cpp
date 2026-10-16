#include <AMReX.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>
#include <iostream>

#include "CRESCENT_Constants.H"
#include "CRESCENT_JonswapSpectrum.H"

namespace {

using amrex::GpuComplex;
using amrex::Real;

bool nearlyEqual(const Real lhs, const Real rhs, const Real tolerance = 2.0e-13)
{
    return std::abs(lhs - rhs) <=
           tolerance * std::max({Real(1.0), std::abs(lhs), std::abs(rhs)});
}

Real magnitude(const GpuComplex<Real>& value)
{
    return std::hypot(value.real(), value.imag());
}

bool isZero(const crescent::init::jonswap::TravelingMode& mode)
{
    return magnitude(mode.eta) == 0.0 && magnitude(mode.phi) == 0.0;
}

bool hasLinearTravelingRelation(
    const crescent::init::jonswap::TravelingMode& mode,
    const Real inverse_gravity,
    const Real omega)
{
    const Real scale = Real(1.0) / (inverse_gravity * omega);
    return nearlyEqual(mode.phi.real(), scale * mode.eta.imag()) &&
           nearlyEqual(mode.phi.imag(), -scale * mode.eta.real());
}

bool hasOppositeTravelingRelation(
    const crescent::init::jonswap::TravelingMode& mode,
    const Real inverse_gravity,
    const Real omega)
{
    const Real scale = Real(1.0) / (inverse_gravity * omega);
    return nearlyEqual(mode.phi.real(), -scale * mode.eta.imag()) &&
           nearlyEqual(mode.phi.imag(), scale * mode.eta.real());
}

int runTest()
{
    namespace jonswap = crescent::init::jonswap;
    using crescent::constants::PI;

    jonswap::SpectrumParameters parameters;
    parameters.nx = 64;
    parameters.dk_x = 0.1;
    parameters.dk_y = 0.1;
    parameters.inverse_gravity = Real(1.0) / Real(9.81);
    parameters.alpha = 0.0081;
    parameters.peak_frequency = 1.0;
    parameters.initial_angle = 0.25;
    parameters.angle_spread = 0.5;
    parameters.directional_exponent = 4;
    parameters.directional_normalization =
        jonswap::cosinePowerNormalization(
            parameters.angle_spread,
            parameters.directional_exponent);
    parameters.amplitude_scale = 1.25;
    parameters.use_directional_spectrum = 1;

    bool ok = true;

    // The general cosine-power normalization integrates to unity and reduces
    // to the historical 2/spread expression only for exponent two.
    const Real quadratic_normalization =
        jonswap::cosinePowerNormalization(parameters.angle_spread, 2);
    const Real normalized_integral =
        parameters.directional_normalization * parameters.angle_spread / PI *
        jonswap::cosinePowerIntegral(parameters.directional_exponent);
    if (!nearlyEqual(
            quadratic_normalization,
            Real(2.0) / parameters.angle_spread) ||
        !nearlyEqual(normalized_integral, 1.0)) {
        std::cerr << "cosine-power spreading is not normalized\n";
        ok = false;
    }
    if (!std::isfinite(jonswap::cosinePowerNormalization(0.5, 1000))) {
        std::cerr << "large directional exponent overflowed normalization\n";
        ok = false;
    }

    auto support_parameters = parameters;
    support_parameters.initial_angle = 0.0;
    support_parameters.angle_spread = 0.5;
    support_parameters.directional_normalization =
        jonswap::cosinePowerNormalization(
            support_parameters.angle_spread,
            support_parameters.directional_exponent);
    if (!jonswap::hasDirectionalSupport(0.8, 0.1, support_parameters) ||
        jonswap::hasDirectionalSupport(-0.8, -0.1, support_parameters)) {
        std::cerr << "directional-support detection is incorrect\n";
        ok = false;
    }

    // A narrow distribution centered at a positive signed angle must not
    // mirror its density onto the negative-ky mode.
    const auto signed_modes = jonswap::makeModePair(
        8, 2, parameters, 0.3, 1.1);
    if (!(magnitude(signed_modes.positive_ky.eta) > 0.0) ||
        !isZero(signed_modes.negative_ky)) {
        std::cerr << "signed directional energy was mirrored across ky=0\n";
        ok = false;
    }

    const Real signed_wavenumber = std::hypot(
        Real(8) * parameters.dk_x,
        Real(2) * parameters.dk_y);
    const Real signed_omega =
        std::sqrt(signed_wavenumber / parameters.inverse_gravity);
    if (!hasLinearTravelingRelation(
            signed_modes.positive_ky,
            parameters.inverse_gravity,
            signed_omega)) {
        std::cerr << "directional mode violates phi_k=-i(g/omega)eta_k\n";
        ok = false;
    }

    // An R2C coefficient represents the conjugate wavevector too, so a
    // distribution centered in the backward half-plane is stored as an
    // opposite-traveling contribution with the required +i(g/omega) relation
    // at the retained +kx coefficient.
    auto backward_parameters = parameters;
    backward_parameters.initial_angle = Real(3.0) * PI / Real(4.0);
    backward_parameters.angle_spread = 0.5;
    backward_parameters.directional_exponent = 2;
    backward_parameters.directional_normalization =
        jonswap::cosinePowerNormalization(
            backward_parameters.angle_spread,
            backward_parameters.directional_exponent);
    const auto represented_backward = jonswap::makeModePair(
        8, 8, backward_parameters, 0.2, 0.4, 0.6, 0.8);
    const Real backward_wavenumber = std::hypot(
        Real(8) * backward_parameters.dk_x,
        Real(8) * backward_parameters.dk_y);
    const Real backward_omega =
        std::sqrt(backward_wavenumber / backward_parameters.inverse_gravity);
    if (!isZero(represented_backward.positive_ky) ||
        !(magnitude(represented_backward.negative_ky.eta) > 0.0) ||
        !hasOppositeTravelingRelation(
            represented_backward.negative_ky,
            backward_parameters.inverse_gravity,
            backward_omega)) {
        std::cerr << "backward directional JONSWAP mode was not represented\n";
        ok = false;
    }

    // When the opposite direction has zero density, the stored coefficient
    // is the forward traveling mode alone, whatever the opposite phase.
    const auto forward_only = jonswap::makeTravelingMode(
        jonswap::normalizedCoefficientMagnitude(
            jonswap::directionalWavenumberDensity(
                Real(8) * parameters.dk_x,
                Real(2) * parameters.dk_y,
                parameters),
            parameters.dk_x * parameters.dk_y,
            parameters.amplitude_scale),
        parameters.inverse_gravity,
        signed_omega,
        0.3);
    const auto full_forward_modes = jonswap::makeModePair(
        8, 2, parameters, 0.3, 1.1, 0.7, 1.9);
    if (!nearlyEqual(
            full_forward_modes.positive_ky.eta.real(),
            forward_only.eta.real()) ||
        !nearlyEqual(
            full_forward_modes.positive_ky.eta.imag(),
            forward_only.eta.imag()) ||
        !nearlyEqual(
            full_forward_modes.positive_ky.phi.real(),
            forward_only.phi.real()) ||
        !nearlyEqual(
            full_forward_modes.positive_ky.phi.imag(),
            forward_only.phi.imag()) ||
        !isZero(full_forward_modes.negative_ky)) {
        std::cerr << "zero opposite density changed the forward spectrum\n";
        ok = false;
    }

    // On the kx=0 real-FFT plane, ky partners must be exact conjugates.
    auto self_plane_parameters = parameters;
    self_plane_parameters.initial_angle = PI / Real(2.0);
    self_plane_parameters.angle_spread = 0.8;
    self_plane_parameters.directional_exponent = 2;
    self_plane_parameters.directional_normalization =
        jonswap::cosinePowerNormalization(
            self_plane_parameters.angle_spread,
            self_plane_parameters.directional_exponent);
    const auto self_plane_modes = jonswap::makeModePair(
        0, 3, self_plane_parameters, 0.7, 2.2);
    if (!(magnitude(self_plane_modes.positive_ky.eta) > 0.0) ||
        !nearlyEqual(
            self_plane_modes.negative_ky.eta.real(),
            self_plane_modes.positive_ky.eta.real()) ||
        !nearlyEqual(
            self_plane_modes.negative_ky.eta.imag(),
            -self_plane_modes.positive_ky.eta.imag()) ||
        !nearlyEqual(
            self_plane_modes.negative_ky.phi.real(),
            self_plane_modes.positive_ky.phi.real()) ||
        !nearlyEqual(
            self_plane_modes.negative_ky.phi.imag(),
            -self_plane_modes.positive_ky.phi.imag())) {
        std::cerr << "kx=0 plane does not satisfy Hermitian ky constraints\n";
        ok = false;
    }

    auto minus_y_parameters = self_plane_parameters;
    minus_y_parameters.initial_angle = -PI / Real(2.0);
    const auto minus_y_modes = jonswap::makeModePair(
        0, 3, minus_y_parameters, 0.7, 2.2, 1.3, 2.7);
    const Real minus_y_wavenumber =
        Real(3) * minus_y_parameters.dk_y;
    const Real minus_y_omega =
        std::sqrt(minus_y_wavenumber / minus_y_parameters.inverse_gravity);
    if (!(magnitude(minus_y_modes.positive_ky.eta) > 0.0) ||
        !hasOppositeTravelingRelation(
            minus_y_modes.positive_ky,
            minus_y_parameters.inverse_gravity,
            minus_y_omega) ||
        !nearlyEqual(
            minus_y_modes.negative_ky.eta.real(),
            minus_y_modes.positive_ky.eta.real()) ||
        !nearlyEqual(
            minus_y_modes.negative_ky.eta.imag(),
            -minus_y_modes.positive_ky.eta.imag())) {
        std::cerr << "negative-y JONSWAP mode violates R2C conjugacy\n";
        ok = false;
    }

    // The even-grid x-Nyquist plane is deliberately empty because it is
    // removed by the initialization filter and is self-conjugate in x.
    const auto nyquist_modes = jonswap::makeModePair(
        parameters.nx / 2, 1, parameters, 0.2, 0.9);
    if (!isZero(nyquist_modes.positive_ky) ||
        !isZero(nyquist_modes.negative_ky)) {
        std::cerr << "x-Nyquist plane was not zeroed\n";
        ok = false;
    }

    // Switching off directional spreading creates a one-dimensional +x
    // spectrum, with no energy at nonzero ky.
    auto long_crested_parameters = parameters;
    long_crested_parameters.use_directional_spectrum = 0;
    const auto transverse_modes = jonswap::makeModePair(
        8, 2, long_crested_parameters, 0.4, 1.2);
    const auto axial_modes = jonswap::makeModePair(
        8, 0, long_crested_parameters, 0.4, 1.2);
    if (!isZero(transverse_modes.positive_ky) ||
        !isZero(transverse_modes.negative_ky) ||
        !(magnitude(axial_modes.positive_ky.eta) > 0.0)) {
        std::cerr << "directional switch-off is not long-crested\n";
        ok = false;
    }

    // Independently reconstruct the expected normalized coefficient from
    // S_k = S_omega g/(2 omega).
    const Real k_x = Real(8) * long_crested_parameters.dk_x;
    const Real omega =
        std::sqrt(k_x / long_crested_parameters.inverse_gravity);
    const Real sigma = omega < long_crested_parameters.peak_frequency
        ? Real(0.07)
        : Real(0.09);
    const Real normalized_offset =
        (omega - long_crested_parameters.peak_frequency) /
        (sigma * long_crested_parameters.peak_frequency);
    const Real gravity = Real(1.0) / long_crested_parameters.inverse_gravity;
    const Real frequency_density =
        long_crested_parameters.alpha * gravity * gravity /
        std::pow(omega, 5) *
        std::exp(-Real(1.25) *
                 std::pow(
                     omega / long_crested_parameters.peak_frequency,
                     -4.0)) *
        std::pow(Real(3.3),
                 std::exp(-Real(0.5) *
                          normalized_offset * normalized_offset));
    const Real wavenumber_density =
        frequency_density * gravity / (Real(2.0) * omega);
    const Real expected_magnitude =
        Real(0.5) * long_crested_parameters.amplitude_scale *
        std::sqrt(
            Real(2.0) * wavenumber_density *
            long_crested_parameters.dk_x);
    if (!nearlyEqual(
            magnitude(axial_modes.positive_ky.eta),
            expected_magnitude) ||
        !hasLinearTravelingRelation(
            axial_modes.positive_ky,
            long_crested_parameters.inverse_gravity,
            omega)) {
        std::cerr << "long-crested coefficient normalization is incorrect\n";
        ok = false;
    }

    return ok ? 0 : 1;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = runTest();
    amrex::Finalize();
    return result;
}

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <iomanip>
#include <iostream>
#include <sstream>

#include "solver/CRESCENT_CheckpointConfigIO.H"

namespace {

template <typename L, typename R>
bool checkEqual(const char* name, const L& lhs, const R& rhs)
{
    if (lhs == rhs) {
        return true;
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "checkpoint config mismatch for " << name
                  << ": " << lhs << " != " << rhs << '\n';
    }
    return false;
}

#define CHECK_FIELD(object_a, object_b, field) \
    ok = checkEqual(#field, (object_a).field, (object_b).field) && ok

int runTest()
{
    crescent::shared::WaveInitializationState initialization;
    initialization.wind_speed_10m = 12.5;
    initialization.fetch = 34567.0;
    initialization.alpha = 0.125;
    initialization.peak_frequency = 0.875;
    initialization.peak_wavelength = 42.5;
    initialization.wave_steepness = 0.1875;
    initialization.initial_angle = -0.25;
    initialization.angle_spread = 0.375;
    initialization.wave_type = "checkpoint test wave";
    initialization.jonswap_random_seed = 987654321;
    initialization.directional_exponent = 6;
    initialization.use_directional_spectrum = 0;
    initialization.jonswap_amplitude_scale = 1.375;
    initialization.wave_mode_x = 9;
    initialization.wave_angle = 0.625;
    initialization.crescent_eps = 0.1875;
    initialization.crescent_kx_ratio = 1.75;
    initialization.crescent_ky_ratio = 1.125;
    initialization.crescent_stokes_order = 4;

    crescent::shared::WaveSpectralState wave;
    wave.perturbation_order = 4;
    wave.beta1_filter = 7.25;
    wave.beta1_filter_l1 = 6.5;
    wave.beta2_filter = 11.75;
    wave.l1_rolloff_alpha = 1234.5;
    wave.l1_rolloff_power = 18;
    wave.adjustment_periods = 2.25;
    wave.peak_wavenumber_override = 0.375;
    wave.peak_wavenumber_physical = 0.3125;

    crescent::nested::AlgorithmConfig nested;
    nested.use_crescent_algorithm = 1;
    nested.l1_box_x_lo = 1.25;
    nested.l1_box_x_hi = 7.75;
    nested.l1_box_y_lo = 2.5;
    nested.l1_box_y_hi = 6.5;
    nested.roi_accurate_fraction = 0.375;
    nested.roi_spectral_fraction = 0.125;
    nested.state_interpolation_order = 1;
    nested.derivative_interpolation_order = 1;
    nested.c3_boundary_penalty = 0.0625;
    nested.c3_coefficient_regularization = 3.5;
    nested.l1_coupling_mode = 1;

    std::stringstream checkpoint;
    checkpoint << std::setprecision(17);
    crescent::solver::checkpoint::config_io::writeWaveInitializationState(
        checkpoint, initialization);
    crescent::solver::checkpoint::config_io::writeWaveSpectralState(
        checkpoint, wave);
    crescent::solver::checkpoint::config_io::writeNestedConfig(
        checkpoint, nested);

    crescent::shared::WaveInitializationState restored_initialization;
    crescent::shared::WaveSpectralState restored_wave;
    crescent::nested::AlgorithmConfig restored_nested;
    crescent::solver::checkpoint::config_io::readWaveInitializationState(
        checkpoint, restored_initialization);
    crescent::solver::checkpoint::config_io::readWaveSpectralState(
        checkpoint, restored_wave);
    crescent::solver::checkpoint::config_io::readNestedConfig(
        checkpoint, restored_nested);

    bool ok = true;
    CHECK_FIELD(initialization, restored_initialization, wind_speed_10m);
    CHECK_FIELD(initialization, restored_initialization, fetch);
    CHECK_FIELD(initialization, restored_initialization, alpha);
    CHECK_FIELD(initialization, restored_initialization, peak_frequency);
    CHECK_FIELD(initialization, restored_initialization, peak_wavelength);
    CHECK_FIELD(initialization, restored_initialization, wave_steepness);
    CHECK_FIELD(initialization, restored_initialization, initial_angle);
    CHECK_FIELD(initialization, restored_initialization, angle_spread);
    CHECK_FIELD(initialization, restored_initialization, wave_type);
    CHECK_FIELD(initialization, restored_initialization, jonswap_random_seed);
    CHECK_FIELD(initialization, restored_initialization, directional_exponent);
    CHECK_FIELD(initialization, restored_initialization, use_directional_spectrum);
    CHECK_FIELD(initialization, restored_initialization, jonswap_amplitude_scale);
    CHECK_FIELD(initialization, restored_initialization, wave_mode_x);
    CHECK_FIELD(initialization, restored_initialization, wave_angle);
    CHECK_FIELD(initialization, restored_initialization, crescent_eps);
    CHECK_FIELD(initialization, restored_initialization, crescent_kx_ratio);
    CHECK_FIELD(initialization, restored_initialization, crescent_ky_ratio);
    CHECK_FIELD(initialization, restored_initialization, crescent_stokes_order);

    CHECK_FIELD(wave, restored_wave, perturbation_order);
    CHECK_FIELD(wave, restored_wave, beta1_filter);
    CHECK_FIELD(wave, restored_wave, beta1_filter_l1);
    CHECK_FIELD(wave, restored_wave, beta2_filter);
    CHECK_FIELD(wave, restored_wave, l1_rolloff_alpha);
    CHECK_FIELD(wave, restored_wave, l1_rolloff_power);
    CHECK_FIELD(wave, restored_wave, adjustment_periods);
    CHECK_FIELD(wave, restored_wave, peak_wavenumber_override);
    CHECK_FIELD(wave, restored_wave, peak_wavenumber_physical);

    CHECK_FIELD(nested, restored_nested, use_crescent_algorithm);
    CHECK_FIELD(nested, restored_nested, l1_box_x_lo);
    CHECK_FIELD(nested, restored_nested, l1_box_x_hi);
    CHECK_FIELD(nested, restored_nested, l1_box_y_lo);
    CHECK_FIELD(nested, restored_nested, l1_box_y_hi);
    CHECK_FIELD(nested, restored_nested, roi_accurate_fraction);
    CHECK_FIELD(nested, restored_nested, roi_spectral_fraction);
    CHECK_FIELD(nested, restored_nested, state_interpolation_order);
    CHECK_FIELD(nested, restored_nested, derivative_interpolation_order);
    CHECK_FIELD(nested, restored_nested, c3_boundary_penalty);
    CHECK_FIELD(nested, restored_nested, c3_coefficient_regularization);
    CHECK_FIELD(nested, restored_nested, l1_coupling_mode);

    return ok ? 0 : 1;
}

#undef CHECK_FIELD

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = runTest();
    amrex::Finalize();
    return result;
}

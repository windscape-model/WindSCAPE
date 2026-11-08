#include "CRESCENT_CheckpointConfigIO.H"

#include <istream>
#include <iomanip>
#include <ostream>
#include <string>

#include <AMReX.H>

namespace {

void requireToken(std::istream& is, const char* expected)
{
    std::string token;
    is >> token;
    if (token != expected) {
        amrex::Abort(
            "Malformed CRESCENT checkpoint header: expected token '" +
            std::string(expected) + "', read '" + token + "'");
    }
}

// Slot of an option retired by the 2026-10 code review. It is written with
// the only supported value, so the format-10 header layout is unchanged, and
// a checkpoint that holds another value cannot be restarted.
template <typename T>
void readRetiredSlot(std::istream& is, const char* name, T supported)
{
    T value{};
    is >> value;
    if (!(value == supported)) {
        amrex::Abort(
            std::string("CRESCENT checkpoint stores a value of the retired "
                        "option crescent.") + name +
            " that this version does not support");
    }
}

} // namespace

namespace crescent::solver::checkpoint::config_io {

void writeWaveInitializationState(
    std::ostream& os,
    const crescent::shared::WaveInitializationState& initialization)
{
    os << "wave_initialization"
       << ' ' << initialization.wind_speed_10m
       << ' ' << initialization.fetch
       << ' ' << initialization.alpha
       << ' ' << initialization.peak_frequency
       << ' ' << initialization.peak_wavelength
       << ' ' << initialization.wave_steepness
       << ' ' << 1  // retired: swell_wavenumber
       << ' ' << initialization.initial_angle
       << ' ' << initialization.angle_spread
       << ' ' << std::quoted(initialization.wave_type)
       << ' ' << initialization.jonswap_random_seed
       << ' ' << initialization.directional_exponent
       << ' ' << initialization.use_directional_spectrum
       << " 1 1"  // retired: correct_jonswap_initialization, full_directional_jonswap
       << ' ' << initialization.jonswap_amplitude_scale
       << ' ' << initialization.wave_mode_x
       << ' ' << initialization.wave_angle
       << ' ' << 1  // retired: correct_stokes_wavevector
       << ' ' << initialization.crescent_eps
       << ' ' << initialization.crescent_kx_ratio
       << ' ' << initialization.crescent_ky_ratio
       << ' ' << initialization.crescent_stokes_order
       << '\n';
}

void readWaveInitializationState(
    std::istream& is,
    crescent::shared::WaveInitializationState& initialization)
{
    requireToken(is, "wave_initialization");
    is >> initialization.wind_speed_10m
       >> initialization.fetch
       >> initialization.alpha
       >> initialization.peak_frequency
       >> initialization.peak_wavelength
       >> initialization.wave_steepness;
    readRetiredSlot(is, "swell_wavenumber", amrex::Real(1.0));
    is >> initialization.initial_angle
       >> initialization.angle_spread
       >> std::quoted(initialization.wave_type)
       >> initialization.jonswap_random_seed
       >> initialization.directional_exponent
       >> initialization.use_directional_spectrum;
    readRetiredSlot(is, "correct_jonswap_initialization", 1);
    readRetiredSlot(is, "full_directional_jonswap", 1);
    is >> initialization.jonswap_amplitude_scale
       >> initialization.wave_mode_x
       >> initialization.wave_angle;
    readRetiredSlot(is, "correct_stokes_wavevector", 1);
    is >> initialization.crescent_eps
       >> initialization.crescent_kx_ratio
       >> initialization.crescent_ky_ratio
       >> initialization.crescent_stokes_order;
}

void writeWaveSpectralState(
    std::ostream& os,
    const crescent::shared::WaveSpectralState& spectral)
{
    os << "wave_spectral"
       << ' ' << spectral.perturbation_order
       << ' ' << spectral.beta1_filter
       << ' ' << spectral.beta1_filter_l1
       << ' ' << spectral.beta2_filter
       << ' ' << spectral.l1_rolloff_alpha
       << ' ' << spectral.l1_rolloff_power
       << ' ' << 1  // retired: l1_square_dealias
       << ' ' << spectral.adjustment_periods
       << ' ' << spectral.peak_wavenumber_override
       << ' ' << 1  // retired: honor_peak_wavenumber_override
       << ' ' << spectral.peak_wavenumber_physical
       << '\n';
}

void readWaveSpectralState(
    std::istream& is,
    crescent::shared::WaveSpectralState& spectral)
{
    requireToken(is, "wave_spectral");
    is >> spectral.perturbation_order
       >> spectral.beta1_filter
       >> spectral.beta1_filter_l1
       >> spectral.beta2_filter
       >> spectral.l1_rolloff_alpha
       >> spectral.l1_rolloff_power;
    readRetiredSlot(is, "l1_square_dealias", 1);
    is >> spectral.adjustment_periods
       >> spectral.peak_wavenumber_override;
    readRetiredSlot(is, "honor_peak_wavenumber_override", 1);
    is >> spectral.peak_wavenumber_physical;
}

void writeNestedConfig(
    std::ostream& os,
    const crescent::nested::AlgorithmConfig& config)
{
    os << "nested_config"
       << ' ' << config.use_crescent_algorithm
       << ' ' << config.l1_box_x_lo
       << ' ' << config.l1_box_x_hi
       << ' ' << config.l1_box_y_lo
       << ' ' << config.l1_box_y_hi
       << ' ' << config.roi_accurate_fraction
       << ' ' << config.roi_spectral_fraction
       << ' ' << 2  // retired: damping_steepness
       << ' ' << 2  // retired: blending_steepness (no effect)
       << ' ' << 1  // retired: blending_tensor_alpha
       << ' ' << 1  // retired: blending_superellipse_gamma
       << ' ' << 4  // retired: blending_superellipse_p
       << ' ' << 0  // retired: blending_correction_y_smoothing_passes
       << ' ' << config.state_interpolation_order
       << ' ' << config.derivative_interpolation_order
       << ' ' << 1  // retired: l1_correct_bilinear_interpolation
       << ' ' << 3  // retired: blending_function_mode
       << ' ' << 0  // retired: use_sequential_blending
       << ' ' << 1  // retired: correct_sequential_blending
       << ' ' << 1   // retired: use_c3_repair
       << ' ' << 1   // retired: c3_apply_x
       << ' ' << 1   // retired: c3_apply_y
       << ' ' << 1   // retired: c3_use_interface
       << ' ' << 24  // retired: c3_buffer_width
       << ' ' << 1   // retired: c3_directional_buffer_widths
       << ' ' << config.c3_boundary_penalty
       << ' ' << config.c3_coefficient_regularization
       << " 0 0 0 0 0 0 0 0 0 0 0 0"  // retired exploratory C3 options
       << ' ' << 1  // retired: filter_periodic_copies_before_linear_derivatives
       << ' ' << 0  // retired: linear_derivative_spectral_filter_mode
       << ' ' << 1  // retired: l1_filter_preserve_mean
       << ' ' << config.l1_coupling_mode
       << ' ' << 0  // retired: l1_rhs_eta_source
       << " 1 1 1 1 1"  // retired CCF consistency switches
       << ' ' << 1  // retired: consistent_applied_pressure
       << ' ' << 1  // retired: l1_final_dno_projection
       << ' ' << 0  // retired: l1_stage_phi_z_cleanup
       << ' ' << 0  // retired: enable_post_timestep_low_pass_filter
       << ' ' << 1  // retired: enable_post_rk4_state_filter_and_blend
       << ' ' << 0  // retired: enable_final_output_blended_derivative_filter
       << ' ' << 1  // retired: consistent_output_diagnostics
       << ' ' << 1  // retired: ccf_band_restrict (the band projection is part of CCF)
       << '\n';
}

void readNestedConfig(
    std::istream& is,
    crescent::nested::AlgorithmConfig& config)
{
    requireToken(is, "nested_config");
    is >> config.use_crescent_algorithm
       >> config.l1_box_x_lo
       >> config.l1_box_x_hi
       >> config.l1_box_y_lo
       >> config.l1_box_y_hi
       >> config.roi_accurate_fraction
       >> config.roi_spectral_fraction;
    readRetiredSlot(is, "damping_steepness", amrex::Real(2.0));
    amrex::Real ignored_blending_steepness = 0.0;  // never had an effect
    is >> ignored_blending_steepness;
    readRetiredSlot(is, "blending_tensor_alpha", amrex::Real(1.0));
    readRetiredSlot(is, "blending_superellipse_gamma", amrex::Real(1.0));
    readRetiredSlot(is, "blending_superellipse_p", amrex::Real(4.0));
    readRetiredSlot(is, "blending_correction_y_smoothing_passes", 0);
    is >> config.state_interpolation_order
       >> config.derivative_interpolation_order;
    readRetiredSlot(is, "l1_correct_bilinear_interpolation", 1);
    readRetiredSlot(is, "blending_function_mode", 3);
    readRetiredSlot(is, "use_sequential_blending", 0);
    readRetiredSlot(is, "correct_sequential_blending", 1);
    readRetiredSlot(is, "use_c3_repair", 1);
    readRetiredSlot(is, "c3_apply_x", 1);
    readRetiredSlot(is, "c3_apply_y", 1);
    readRetiredSlot(is, "c3_use_interface", 1);
    readRetiredSlot(is, "c3_buffer_width", 24);
    readRetiredSlot(is, "c3_directional_buffer_widths", 1);
    is >> config.c3_boundary_penalty
       >> config.c3_coefficient_regularization;
    for (const char* name : {"c3_x_pass_y_smoothing_passes", "c3_x_pass_corner_taper",
                             "c3_x_alpha_y_smoothing_passes"}) {
        readRetiredSlot(is, name, 0);
    }
    for (const char* name : {"c3_x_alpha_y_regularization", "c3_y_alpha_x_regularization",
                             "c3_x_alpha_y_mean_regularization",
                             "c3_y_alpha_x_mean_regularization"}) {
        readRetiredSlot(is, name, amrex::Real(0.0));
    }
    for (const char* name : {"c3_repair_after_filter", "c3_x_transverse_harmonics",
                             "c3_y_transverse_harmonics",
                             "c3_x_alpha_y_projection_harmonics",
                             "c3_y_alpha_x_projection_harmonics"}) {
        readRetiredSlot(is, name, 0);
    }
    readRetiredSlot(is, "filter_periodic_copies_before_linear_derivatives", 1);
    readRetiredSlot(is, "linear_derivative_spectral_filter_mode", 0);
    readRetiredSlot(is, "l1_filter_preserve_mean", 1);
    is >> config.l1_coupling_mode;
    if (config.l1_coupling_mode != 1 && config.l1_coupling_mode != 2) {
        amrex::Abort(
            "CRESCENT checkpoint uses l1_coupling_mode " +
            std::to_string(config.l1_coupling_mode) +
            "; this version supports 1 and 2 (mode 0 was retired)");
    }
    readRetiredSlot(is, "l1_rhs_eta_source", 0);
    for (const char* name : {"ccf_background_tendency", "ccf_linear_consistency",
                             "ccf_adjustment_consistency", "ccf_transform_admissibility",
                             "ccf_w1_projection_consistency"}) {
        readRetiredSlot(is, name, 1);
    }
    readRetiredSlot(is, "consistent_applied_pressure", 1);
    readRetiredSlot(is, "l1_final_dno_projection", 1);
    readRetiredSlot(is, "l1_stage_phi_z_cleanup", 0);
    readRetiredSlot(is, "enable_post_timestep_low_pass_filter", 0);
    readRetiredSlot(is, "enable_post_rk4_state_filter_and_blend", 1);
    readRetiredSlot(is, "enable_final_output_blended_derivative_filter", 0);
    readRetiredSlot(is, "consistent_output_diagnostics", 1);
    readRetiredSlot(is, "ccf_band_restrict", 1);
}

} // namespace crescent::solver::checkpoint::config_io

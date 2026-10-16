#include "CRESCENT_SetupParameters.H"

#include <cmath>
#include <string>

#include <AMReX_ParmParse.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_NestedAlgorithmInputs.H"
#include "CRESCENT_RetiredOptions.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolverState.H"
#include "CRESCENT_WaveState.H"

namespace {

// Run control, output, checkpoint and wave inputs. A key that a deck leaves
// out keeps the default of its state field.
void readRuntimeInputs(
    amrex::ParmParse& pp,
    crescent::shared::RuntimeControlState& runtime,
    crescent::shared::WaveRuntimeState& wave)
{
    amrex::ParmParse top;
    top.query("max_step", runtime.max_step);
    top.query("stop_time", runtime.stop_time);
    top.query("start_time", runtime.start_time);

    auto& reporting = runtime.reporting;
    pp.query("plot_file", reporting.plot_file);
    pp.query("plot_int", reporting.plot_int);
    pp.query("tecplot_int", reporting.tecplot_int);
    pp.query("write_omega_a_edge_diagnostics", reporting.write_omega_a_edge_diagnostics);
    pp.query("omega_a_edge_diagnostic_int", reporting.omega_a_edge_diagnostic_int);
    pp.query("verbosity", runtime.verbosity);
    // NetCDF and HDF5 output were never implemented.
    crescent::shared::requireRetiredValue(pp, "netcdf_int", -1);
    crescent::shared::requireRetiredValue(pp, "hdf5_int", -1);

    auto& checkpoint = runtime.checkpoint;
    pp.query("check_file", checkpoint.file_prefix);
    pp.query("check_int", checkpoint.interval);
    pp.query("restart", checkpoint.restart_file);

    auto& initialization = wave.initialization;
    auto& spectral = wave.spectral;
    pp.query("wind_speed_10m", initialization.wind_speed_10m);
    pp.query("fetch", initialization.fetch);
    pp.query("wave_type", initialization.wave_type);
    pp.query("jonswap_random_seed", initialization.jonswap_random_seed);
    pp.query("directional_exponent", initialization.directional_exponent);
    pp.query("use_directional_spectrum", initialization.use_directional_spectrum);
    pp.query("jonswap_amplitude_scale", initialization.jonswap_amplitude_scale);
    pp.query("wave_steepness", initialization.wave_steepness);
    pp.query("wave_mode_x", initialization.wave_mode_x);
    pp.query("wave_angle", initialization.wave_angle);
    pp.query("crescent_eps", initialization.crescent_eps);
    pp.query("crescent_kx_ratio", initialization.crescent_kx_ratio);
    pp.query("crescent_ky_ratio", initialization.crescent_ky_ratio);
    pp.query("crescent_stokes_order", initialization.crescent_stokes_order);
    pp.query("initial_angle", initialization.initial_angle);
    pp.query("angle_spread", initialization.angle_spread);
    pp.query("beta1_filter", spectral.beta1_filter);
    pp.query("beta1_filter_l1", spectral.beta1_filter_l1);
    pp.query("l1_rolloff_alpha", spectral.l1_rolloff_alpha);
    pp.query("l1_rolloff_power", spectral.l1_rolloff_power);
    pp.query("adjustment_periods", spectral.adjustment_periods);
    pp.query("beta2_filter", spectral.beta2_filter);
    pp.query("peak_wavenumber_override", spectral.peak_wavenumber_override);
    pp.query("perturbation_order", spectral.perturbation_order);
    for (const char* name : {"correct_jonswap_initialization",
                             "full_directional_jonswap",
                             "correct_stokes_wavevector",
                             "l1_square_dealias",
                             "honor_peak_wavenumber_override"}) {
        crescent::shared::requireRetiredValue(pp, name, 1);
    }

    for (const char* name : {"alpha", "peak_frequency", "peak_wavelength"}) {
        if (pp.contains(name)) {
            amrex::Abort(
                std::string("crescent.") + name + " is derived from "
                "crescent.wind_speed_10m and crescent.fetch; it is not a "
                "supported input");
        }
    }
    if (pp.contains("swell_wavenumber")) {
        amrex::Abort(
            "crescent.swell_wavenumber has no implemented wave initializer; "
            "remove this unsupported input");
    }

    // JONSWAP takes its step from the peak period (T_peak/20) unless
    // crescent.dt_override is set; every other initializer needs crescent.dt.
    const bool is_jonswap = initialization.wave_type == "jonswap";
    const bool has_dt = pp.contains("dt");
    const bool has_dt_override = pp.query("dt_override", runtime.dt_override);
    pp.query("dt", runtime.dt[0]);
    if (is_jonswap && has_dt) {
        amrex::Abort(
            "crescent.dt is not used by JONSWAP initialization; use "
            "crescent.dt_override or omit both for T_peak/20");
    }
    if (is_jonswap && has_dt_override &&
        (!(runtime.dt_override > 0.0) ||
         !std::isfinite(runtime.dt_override))) {
        amrex::Abort(
            "an explicit crescent.dt_override for JONSWAP initialization "
            "must be a positive finite number; omit it to use T_peak/20");
    }
    if (!is_jonswap && has_dt_override) {
        amrex::Abort(
            "crescent.dt_override is only used by JONSWAP initialization; "
            "use crescent.dt for this wave_type");
    }
}

void validateOutputInterval(int interval, const char* name)
{
    if (interval < -1) {
        amrex::Abort(std::string(name) + " must be >= -1");
    }
}

void validateRuntimeStepState(const crescent::shared::RuntimeControlState& runtime)
{
    if (runtime.max_step < 0) {
        amrex::Abort("max_step must be >= 0");
    }
    if (!std::isfinite(runtime.start_time)) {
        amrex::Abort("start_time must be finite");
    }
    if (!std::isfinite(runtime.stop_time)) {
        amrex::Abort("stop_time must be finite");
    }
    if (runtime.stop_time < runtime.start_time) {
        amrex::Abort("stop_time must be >= start_time");
    }
    if (!std::isfinite(runtime.dt_override)) {
        amrex::Abort("crescent.dt_override must be finite");
    }
}

void validateReportingState(const crescent::shared::RuntimeReportingState& reporting)
{
    validateOutputInterval(reporting.plot_int, "crescent.plot_int");
    validateOutputInterval(reporting.tecplot_int, "crescent.tecplot_int");

    if ((reporting.plot_int > 0 || reporting.tecplot_int > 0) &&
        reporting.plot_file.empty()) {
        amrex::Abort("crescent.plot_file must be non-empty when plot or Tecplot output is enabled");
    }
    if (reporting.write_omega_a_edge_diagnostics != 0 &&
        reporting.write_omega_a_edge_diagnostics != 1) {
        amrex::Abort("crescent.write_omega_a_edge_diagnostics must be 0 or 1");
    }
    if (reporting.omega_a_edge_diagnostic_int <= 0) {
        amrex::Abort("crescent.omega_a_edge_diagnostic_int must be >= 1");
    }
}

void validateCheckpointState(const crescent::shared::RuntimeCheckpointState& checkpoint)
{
    if (checkpoint.interval != -1 && checkpoint.interval <= 0) {
        amrex::Abort("crescent.check_int must be -1 or >= 1");
    }
    if (checkpoint.interval > 0 && checkpoint.file_prefix.empty()) {
        amrex::Abort("crescent.check_file must be non-empty when checkpointing is enabled");
    }
}

void validateWaveState(const crescent::shared::WaveRuntimeState& wave)
{
    const auto& initialization = wave.initialization;
    const auto& spectral = wave.spectral;

    if (spectral.perturbation_order < 1) {
        amrex::Abort("crescent.perturbation_order must be >= 1");
    }
    if (!(spectral.beta1_filter > 0.0) ||
        !std::isfinite(spectral.beta1_filter)) {
        amrex::Abort("crescent.beta1_filter must be a positive finite number");
    }
    if (!(spectral.beta2_filter > 0.0) ||
        !std::isfinite(spectral.beta2_filter)) {
        amrex::Abort("crescent.beta2_filter must be a positive finite number");
    }
    if (spectral.beta1_filter_l1 != -1.0 &&
        (!(spectral.beta1_filter_l1 > 0.0) ||
         !std::isfinite(spectral.beta1_filter_l1))) {
        amrex::Abort(
            "crescent.beta1_filter_l1 must be -1 or a positive finite number");
    }
    if (!(spectral.l1_rolloff_alpha > 0.0) ||
        !std::isfinite(spectral.l1_rolloff_alpha)) {
        amrex::Abort(
            "crescent.l1_rolloff_alpha must be a positive finite number");
    }
    if (spectral.l1_rolloff_power < 1) {
        amrex::Abort("crescent.l1_rolloff_power must be >= 1");
    }
    if (spectral.adjustment_periods < 0.0 ||
        !std::isfinite(spectral.adjustment_periods)) {
        amrex::Abort(
            "crescent.adjustment_periods must be a non-negative finite number");
    }
    if (spectral.peak_wavenumber_override != -1.0 &&
        (!(spectral.peak_wavenumber_override > 0.0) ||
         !std::isfinite(spectral.peak_wavenumber_override))) {
        amrex::Abort(
            "crescent.peak_wavenumber_override must be -1 or a positive finite number");
    }
    if (initialization.jonswap_random_seed < -1) {
        amrex::Abort("crescent.jonswap_random_seed must be >= -1");
    }
    if (initialization.use_directional_spectrum != 0 &&
        initialization.use_directional_spectrum != 1) {
        amrex::Abort("crescent.use_directional_spectrum must be 0 or 1");
    }
    if (initialization.directional_exponent < 0) {
        amrex::Abort("crescent.directional_exponent must be >= 0");
    }
    if (!(initialization.jonswap_amplitude_scale > 0.0) ||
        !std::isfinite(initialization.jonswap_amplitude_scale)) {
        amrex::Abort("crescent.jonswap_amplitude_scale must be a positive finite number");
    }
    if (initialization.wave_steepness < 0.0 ||
        !std::isfinite(initialization.wave_steepness)) {
        amrex::Abort(
            "crescent.wave_steepness must be a non-negative finite number");
    }
    if (initialization.wave_type == "jonswap") {
        if (initialization.wind_speed_10m <= 0.0 ||
            !std::isfinite(initialization.wind_speed_10m)) {
            amrex::Abort(
                "crescent.wind_speed_10m must be positive and finite for JONSWAP initialization");
        }
        if (initialization.fetch <= 0.0 || !std::isfinite(initialization.fetch)) {
            amrex::Abort(
                "crescent.fetch must be positive and finite for JONSWAP initialization");
        }
        if (initialization.angle_spread < 0.0 ||
            !std::isfinite(initialization.angle_spread)) {
            amrex::Abort(
                "crescent.angle_spread must be a non-negative finite number");
        }
        if (initialization.use_directional_spectrum == 1 &&
            (!(initialization.angle_spread > 0.0) ||
             initialization.angle_spread > crescent::constants::PI2)) {
            amrex::Abort(
                "crescent.angle_spread must be in (0, 2*pi] for directional "
                "JONSWAP initialization");
        }
        if (!std::isfinite(initialization.initial_angle)) {
            amrex::Abort("crescent.initial_angle must be finite");
        }
    } else if (initialization.wave_type == "stokes") {
        if (initialization.wave_mode_x <= 0) {
            amrex::Abort("crescent.wave_mode_x must be >= 1 for Stokes initialization");
        }
        if (!std::isfinite(initialization.wave_angle)) {
            amrex::Abort("crescent.wave_angle must be finite for Stokes initialization");
        }
    } else if (initialization.wave_type == "crescent") {
        if (initialization.wave_mode_x <= 0) {
            amrex::Abort("crescent.wave_mode_x must be >= 1 for crescent initialization");
        }
        if (initialization.crescent_stokes_order < 1 ||
            initialization.crescent_stokes_order > 5) {
            amrex::Abort("crescent.crescent_stokes_order must be in [1, 5]");
        }
        if (!(initialization.crescent_kx_ratio > 0.0) ||
            !std::isfinite(initialization.crescent_kx_ratio) ||
            !(initialization.crescent_ky_ratio > 0.0) ||
            !std::isfinite(initialization.crescent_ky_ratio)) {
            amrex::Abort(
                "crescent.crescent_kx_ratio and crescent_ky_ratio must be positive and finite");
        }
        if (!std::isfinite(initialization.crescent_eps)) {
            amrex::Abort("crescent.crescent_eps must be finite");
        }
    } else {
        amrex::Abort(
            "crescent.wave_type must be 'stokes', 'crescent', or 'jonswap'");
    }
}

void validateRuntimeControlState(
    const crescent::shared::RuntimeControlState& runtime,
    const crescent::shared::WaveRuntimeState& wave)
{
    validateRuntimeStepState(runtime);
    validateReportingState(runtime.reporting);
    validateCheckpointState(runtime.checkpoint);

    // JONSWAP selects its step after the spectrum has supplied the peak
    // period. Every other initializer needs a usable step before setup exits;
    // otherwise the 1e100 allocation sentinel (or a user-supplied zero or
    // negative value) reaches the run loop.
    if (wave.initialization.wave_type != "jonswap" &&
        (runtime.dt.empty() ||
         !(runtime.dt[0] > 0.0) ||
         !std::isfinite(runtime.dt[0]))) {
        amrex::Abort(
            "non-JONSWAP runs require crescent.dt to be a positive finite number");
    }
}

} // namespace

namespace crescent::solver::setup {

void readParameters(
    crescent::SolverState& state,
    int max_level,
    const std::string& solver_pp_prefix)
{
    auto& runtime = state.core.runtime;
    auto& wave = state.core.wave;
    auto& nested_config = state.nested.config;

    amrex::ParmParse pp(solver_pp_prefix);
    readRuntimeInputs(pp, runtime, wave);
    const auto nested_inputs = crescent::nested::config::readAlgorithmInputs(pp);
    nested_config = nested_inputs.config;

    crescent::nested::config::validateAlgorithmSelection(nested_config, max_level);
    crescent::nested::config::validateRefinementRegionInputs(
        nested_config, nested_inputs.user_specified_options, max_level);
    crescent::nested::config::validateInterpolationInputs(nested_config);
    validateWaveState(wave);
    crescent::nested::config::validateBlendingInputs(nested_config);
    crescent::nested::config::validateC3Inputs(nested_config);
    crescent::nested::config::emitInputWarnings(
        nested_config, nested_inputs.user_specified_options);
    validateRuntimeControlState(runtime, wave);
}

} // namespace crescent::solver::setup

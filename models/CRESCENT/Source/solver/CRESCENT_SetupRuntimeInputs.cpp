#include "CRESCENT_SetupRuntimeInputs.H"

#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_WaveState.H"

#include <cmath>
#include <optional>
#include <string>

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

namespace {

using ParsedRuntimeInputs = crescent::solver::setup::runtime_inputs::ParsedRuntimeInputs;
using RuntimeStepInputs = crescent::solver::setup::runtime_inputs::RuntimeStepInputs;
using OutputInputs = crescent::solver::setup::runtime_inputs::OutputInputs;
using CheckpointInputs = crescent::solver::setup::runtime_inputs::CheckpointInputs;
using WaveInputs = crescent::solver::setup::runtime_inputs::WaveInputs;

template <typename T>
std::optional<T> queryOptional(amrex::ParmParse& pp, const char* name)
{
    T value{};
    if (pp.query(name, value)) {
        return value;
    }
    return std::nullopt;
}

template <typename T>
void applyIfSpecified(const std::optional<T>& input, T& target)
{
    if (input.has_value()) {
        target = *input;
    }
}

RuntimeStepInputs readRuntimeParameters()
{
    RuntimeStepInputs inputs;
    amrex::ParmParse pp;
    inputs.max_step = queryOptional<int>(pp, "max_step");
    inputs.stop_time = queryOptional<amrex::Real>(pp, "stop_time");
    inputs.start_time = queryOptional<amrex::Real>(pp, "start_time");
    return inputs;
}

OutputInputs readOutputParameters(amrex::ParmParse& pp)
{
    OutputInputs inputs;
    inputs.plot_file =
        queryOptional<std::string>(pp, "plot_file");
    inputs.plot_int =
        queryOptional<int>(pp, "plot_int");
    inputs.tecplot_int =
        queryOptional<int>(pp, "tecplot_int");
    inputs.netcdf_int =
        queryOptional<int>(pp, "netcdf_int");
    inputs.hdf5_int =
        queryOptional<int>(pp, "hdf5_int");
    inputs.write_omega_a_edge_diagnostics =
        queryOptional<int>(pp, "write_omega_a_edge_diagnostics");
    inputs.omega_a_edge_diagnostic_int =
        queryOptional<int>(pp, "omega_a_edge_diagnostic_int");
    inputs.verbosity =
        queryOptional<int>(pp, "verbosity");
    return inputs;
}

CheckpointInputs readCheckpointParameters(amrex::ParmParse& pp)
{
    CheckpointInputs inputs;
    inputs.check_file =
        queryOptional<std::string>(pp, "check_file");
    inputs.check_int =
        queryOptional<int>(pp, "check_int");
    inputs.restart =
        queryOptional<std::string>(pp, "restart");
    return inputs;
}

WaveInputs readWaveParameters(amrex::ParmParse& pp)
{
    WaveInputs inputs;
    inputs.wind_speed_10m =
        queryOptional<amrex::Real>(pp, "wind_speed_10m");
    inputs.fetch =
        queryOptional<amrex::Real>(pp, "fetch");
    inputs.wave_type =
        queryOptional<std::string>(pp, "wave_type");
    inputs.jonswap_random_seed =
        queryOptional<amrex::Long>(pp, "jonswap_random_seed");

    inputs.directional_exponent =
        queryOptional<int>(pp, "directional_exponent");
    inputs.use_directional_spectrum =
        queryOptional<int>(pp, "use_directional_spectrum");
    inputs.alpha =
        queryOptional<amrex::Real>(pp, "alpha");
    inputs.peak_frequency =
        queryOptional<amrex::Real>(pp, "peak_frequency");
    inputs.peak_wavelength =
        queryOptional<amrex::Real>(pp, "peak_wavelength");
    inputs.wave_steepness =
        queryOptional<amrex::Real>(pp, "wave_steepness");
    inputs.swell_wavenumber =
        queryOptional<amrex::Real>(pp, "swell_wavenumber");
    inputs.wave_mode_x =
        queryOptional<int>(pp, "wave_mode_x");
    inputs.wave_angle =
        queryOptional<amrex::Real>(pp, "wave_angle");
    inputs.initial_angle =
        queryOptional<amrex::Real>(pp, "initial_angle");
    inputs.angle_spread =
        queryOptional<amrex::Real>(pp, "angle_spread");
    inputs.beta1_filter =
        queryOptional<amrex::Real>(pp, "beta1_filter");
    inputs.beta2_filter =
        queryOptional<amrex::Real>(pp, "beta2_filter");
    inputs.peak_wavenumber_override =
        queryOptional<amrex::Real>(pp, "peak_wavenumber_override");
    inputs.perturbation_order =
        queryOptional<int>(pp, "perturbation_order");
    return inputs;
}

void applyRuntimeStepInputs(
    crescent::shared::RuntimeControlState& runtime,
    int finest_level,
    const RuntimeStepInputs& inputs)
{
    applyIfSpecified(inputs.max_step, runtime.max_step);
    applyIfSpecified(inputs.stop_time, runtime.stop_time);
    applyIfSpecified(inputs.start_time, runtime.start_time);
    applyIfSpecified(inputs.dt_override, runtime.dt_override);
    applyIfSpecified(inputs.dt, runtime.dt[0]);

    for (int lev = 1; lev <= finest_level; ++lev) {
        runtime.dt[lev] = runtime.dt[0];
    }
}

void applyOutputInputs(
    crescent::shared::RuntimeControlState& runtime,
    const OutputInputs& inputs)
{
    auto& reporting = runtime.reporting;
    applyIfSpecified(inputs.plot_file, reporting.plot_file);
    applyIfSpecified(inputs.plot_int, reporting.plot_int);
    applyIfSpecified(inputs.tecplot_int, reporting.tecplot_int);
    applyIfSpecified(inputs.netcdf_int, reporting.netcdf_int);
    applyIfSpecified(inputs.hdf5_int, reporting.hdf5_int);
    applyIfSpecified(
        inputs.write_omega_a_edge_diagnostics,
        reporting.write_omega_a_edge_diagnostics);
    applyIfSpecified(
        inputs.omega_a_edge_diagnostic_int,
        reporting.omega_a_edge_diagnostic_int);
    applyIfSpecified(inputs.verbosity, runtime.verbosity);
}

void applyCheckpointInputs(
    crescent::shared::RuntimeControlState& runtime,
    const CheckpointInputs& inputs)
{
    auto& checkpoint = runtime.checkpoint;
    applyIfSpecified(inputs.check_file, checkpoint.file_prefix);
    applyIfSpecified(inputs.check_int, checkpoint.interval);
    applyIfSpecified(inputs.restart, checkpoint.restart_file);
}

void applyWaveInputs(
    crescent::shared::WaveRuntimeState& wave,
    const WaveInputs& inputs)
{
    auto& initialization = wave.initialization;
    auto& spectral = wave.spectral;

    applyIfSpecified(inputs.wind_speed_10m, initialization.wind_speed_10m);
    applyIfSpecified(inputs.fetch, initialization.fetch);
    applyIfSpecified(inputs.wave_type, initialization.wave_type);
    applyIfSpecified(inputs.jonswap_random_seed, initialization.jonswap_random_seed);
    applyIfSpecified(inputs.directional_exponent, initialization.directional_exponent);
    applyIfSpecified(inputs.use_directional_spectrum, initialization.use_directional_spectrum);
    applyIfSpecified(inputs.alpha, initialization.alpha);
    applyIfSpecified(inputs.peak_frequency, initialization.peak_frequency);
    applyIfSpecified(inputs.peak_wavelength, initialization.peak_wavelength);
    applyIfSpecified(inputs.wave_steepness, initialization.wave_steepness);
    applyIfSpecified(inputs.swell_wavenumber, initialization.swell_wavenumber);
    applyIfSpecified(inputs.wave_mode_x, initialization.wave_mode_x);
    applyIfSpecified(inputs.wave_angle, initialization.wave_angle);
    applyIfSpecified(inputs.initial_angle, initialization.initial_angle);
    applyIfSpecified(inputs.angle_spread, initialization.angle_spread);
    applyIfSpecified(inputs.beta1_filter, spectral.beta1_filter);
    applyIfSpecified(inputs.beta2_filter, spectral.beta2_filter);
    applyIfSpecified(inputs.peak_wavenumber_override, spectral.peak_wavenumber_override);
    applyIfSpecified(inputs.perturbation_order, spectral.perturbation_order);

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
    validateOutputInterval(reporting.netcdf_int, "crescent.netcdf_int");
    validateOutputInterval(reporting.hdf5_int, "crescent.hdf5_int");

    if ((reporting.plot_int > 0 || reporting.tecplot_int > 0) &&
        reporting.plot_file.empty()) {
        amrex::Abort("crescent.plot_file must be non-empty when plot or Tecplot output is enabled");
    }
    if (reporting.netcdf_int > 0) {
        amrex::Abort("crescent.netcdf_int requested NetCDF output, but NetCDF output is not implemented");
    }
    if (reporting.hdf5_int > 0) {
        amrex::Abort("crescent.hdf5_int requested HDF5 output, but HDF5 output is not implemented");
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

} // namespace

namespace crescent::solver::setup::runtime_inputs {

ParsedRuntimeInputs readRuntimeInputs(amrex::ParmParse& pp)
{
    ParsedRuntimeInputs inputs;
    inputs.runtime = readRuntimeParameters();
    inputs.output = readOutputParameters(pp);
    inputs.checkpoint = readCheckpointParameters(pp);
    inputs.wave = readWaveParameters(pp);
    inputs.runtime.dt =
        queryOptional<amrex::Real>(pp, "dt");
    inputs.runtime.dt_override =
        queryOptional<amrex::Real>(pp, "dt_override");
    return inputs;
}

void validateWaveState(const crescent::shared::WaveRuntimeState& wave)
{
    const auto& initialization = wave.initialization;
    const auto& spectral = wave.spectral;

    if (spectral.perturbation_order < 1) {
        amrex::Abort("crescent.perturbation_order must be >= 1");
    }
    if (initialization.jonswap_random_seed < -1) {
        amrex::Abort("crescent.jonswap_random_seed must be >= -1");
    }
    if (initialization.wave_steepness < 0.0) {
        amrex::Abort("crescent.wave_steepness must be >= 0");
    }
    if (initialization.wave_type == "jonswap") {
        if (initialization.wind_speed_10m <= 0.0 ||
            !std::isfinite(initialization.wind_speed_10m)) {
            amrex::Abort("crescent.wind_speed_10m must be positive for JONSWAP initialization");
        }
        if (initialization.fetch <= 0.0 || !std::isfinite(initialization.fetch)) {
            amrex::Abort("crescent.fetch must be positive for JONSWAP initialization");
        }
    } else if (initialization.wave_type == "stokes") {
        if (initialization.wave_mode_x <= 0) {
            amrex::Abort("crescent.wave_mode_x must be >= 1 for Stokes initialization");
        }
    }
}

void validateRuntimeControlState(const crescent::shared::RuntimeControlState& runtime)
{
    validateRuntimeStepState(runtime);
    validateReportingState(runtime.reporting);
    validateCheckpointState(runtime.checkpoint);
}

void applyRuntimeInputs(
    crescent::shared::RuntimeControlState& runtime,
    crescent::shared::WaveRuntimeState& wave,
    int finest_level,
    const ParsedRuntimeInputs& inputs)
{
    applyRuntimeStepInputs(runtime, finest_level, inputs.runtime);
    applyOutputInputs(runtime, inputs.output);
    applyCheckpointInputs(runtime, inputs.checkpoint);
    applyWaveInputs(wave, inputs.wave);
}

} // namespace crescent::solver::setup::runtime_inputs

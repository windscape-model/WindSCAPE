#pragma once

#include <optional>
#include <string>

#include <AMReX_INT.H>
#include <AMReX_REAL.H>

namespace amrex {

class ParmParse;

} // namespace amrex

namespace crescent::shared {

struct RuntimeControlState;
struct WaveRuntimeState;

} // namespace crescent::shared

namespace crescent::solver::setup::runtime_inputs {

struct RuntimeStepInputs
{
    std::optional<int> max_step;
    std::optional<amrex::Real> start_time;
    std::optional<amrex::Real> stop_time;
    std::optional<amrex::Real> dt;
    std::optional<amrex::Real> dt_override;
};

struct OutputInputs
{
    std::optional<int> verbosity;
    std::optional<std::string> plot_file;
    std::optional<int> plot_int;
    std::optional<int> tecplot_int;
    std::optional<int> netcdf_int;
    std::optional<int> hdf5_int;
    std::optional<int> write_omega_a_edge_diagnostics;
    std::optional<int> omega_a_edge_diagnostic_int;
};

struct CheckpointInputs
{
    std::optional<std::string> check_file;
    std::optional<int> check_int;
    std::optional<std::string> restart;
};

struct WaveInputs
{
    std::optional<int> perturbation_order;
    std::optional<amrex::Real> wind_speed_10m;
    std::optional<amrex::Real> fetch;
    std::optional<amrex::Real> alpha;
    std::optional<amrex::Real> peak_frequency;
    std::optional<amrex::Real> peak_wavelength;
    std::optional<amrex::Real> wave_steepness;
    std::optional<amrex::Real> swell_wavenumber;
    std::optional<amrex::Real> initial_angle;
    std::optional<amrex::Real> angle_spread;
    std::optional<amrex::Real> beta1_filter;
    std::optional<amrex::Real> beta2_filter;
    std::optional<amrex::Real> peak_wavenumber_override;
    std::optional<std::string> wave_type;
    std::optional<amrex::Long> jonswap_random_seed;
    std::optional<int> directional_exponent;
    std::optional<int> use_directional_spectrum;
    std::optional<int> wave_mode_x;
    std::optional<amrex::Real> wave_angle;
};

struct ParsedRuntimeInputs
{
    RuntimeStepInputs runtime;
    OutputInputs output;
    CheckpointInputs checkpoint;
    WaveInputs wave;
};

ParsedRuntimeInputs readRuntimeInputs(
    amrex::ParmParse& pp);

void validateWaveState(const shared::WaveRuntimeState& wave);

void validateRuntimeControlState(const shared::RuntimeControlState& runtime);

void applyRuntimeInputs(
    shared::RuntimeControlState& runtime,
    shared::WaveRuntimeState& wave,
    int finest_level,
    const ParsedRuntimeInputs& inputs);

} // namespace crescent::solver::setup::runtime_inputs

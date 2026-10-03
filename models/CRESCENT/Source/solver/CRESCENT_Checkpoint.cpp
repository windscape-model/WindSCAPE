#include "CRESCENT_Checkpoint.H"

#include <fstream>
#include <sstream>
#include <string>

#include <AMReX_FileSystem.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_RestartRebuild.H"
#include "CRESCENT_RuntimeReportingContext.H"
#include "CRESCENT_SetupInitialization.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverInitializationContext.H"
#include "CRESCENT_SolverOwnership.H"

namespace {

constexpr int CheckpointFormatVersion = 1;
constexpr int CheckpointFileNameDigits = 5;

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

template <typename T>
void writeVector(std::ostream& os, const char* name, const amrex::Vector<T>& values)
{
    os << name;
    for (const auto& value : values) {
        os << ' ' << value;
    }
    os << '\n';
}

template <typename T>
void readVector(std::istream& is, const char* name, amrex::Vector<T>& values)
{
    requireToken(is, name);
    for (auto& value : values) {
        is >> value;
    }
}

std::string checkpointName(const crescent::shared::RuntimeControlState& runtime)
{
    return amrex::Concatenate(
        runtime.checkpoint.file_prefix,
        runtime.istep[0],
        CheckpointFileNameDigits);
}

void writeGeometryMetadata(
    std::ostream& os,
    int lev,
    const amrex::Geometry& geom)
{
    const auto& domain = geom.Domain();
    os << "level " << lev << '\n';
    os << "geometry_domain";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << domain.smallEnd(d);
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << domain.bigEnd(d);
    }
    os << '\n';

    os << "geometry_prob_lo";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << geom.ProbLo(d);
    }
    os << '\n';

    os << "geometry_prob_hi";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << geom.ProbHi(d);
    }
    os << '\n';

    os << "geometry_periodic";
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        os << ' ' << geom.isPeriodic(d);
    }
    os << '\n';
}

void validateGeometryMetadata(
    std::istream& is,
    int expected_lev,
    const amrex::Geometry& geom)
{
    requireToken(is, "level");
    int lev = -1;
    is >> lev;
    if (lev != expected_lev) {
        amrex::Abort("CRESCENT checkpoint level metadata is out of order");
    }

    requireToken(is, "geometry_domain");
    amrex::IntVect small;
    amrex::IntVect big;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        is >> small[d];
    }
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        is >> big[d];
    }
    if (amrex::Box(small, big) != geom.Domain()) {
        amrex::Abort("CRESCENT restart geometry domain does not match the checkpoint");
    }

    requireToken(is, "geometry_prob_lo");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        amrex::Real value = 0.0;
        is >> value;
        if (value != geom.ProbLo(d)) {
            amrex::Abort("CRESCENT restart geometry ProbLo does not match the checkpoint");
        }
    }

    requireToken(is, "geometry_prob_hi");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        amrex::Real value = 0.0;
        is >> value;
        if (value != geom.ProbHi(d)) {
            amrex::Abort("CRESCENT restart geometry ProbHi does not match the checkpoint");
        }
    }

    requireToken(is, "geometry_periodic");
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        int value = 0;
        is >> value;
        if (value != static_cast<int>(geom.isPeriodic(d))) {
            amrex::Abort("CRESCENT restart periodicity does not match the checkpoint");
        }
    }
}

void writeNestedConfig(std::ostream& os, const crescent::nested::AlgorithmConfig& config)
{
    os << "nested_config"
       << ' ' << config.use_crescent_algorithm
       << ' ' << config.l1_box_x_lo
       << ' ' << config.l1_box_x_hi
       << ' ' << config.l1_box_y_lo
       << ' ' << config.l1_box_y_hi
       << ' ' << config.roi_accurate_fraction
       << ' ' << config.roi_spectral_fraction
       << ' ' << config.damping_steepness
       << ' ' << config.blending_steepness
       << ' ' << config.blending_tensor_alpha
       << ' ' << config.blending_superellipse_gamma
       << ' ' << config.blending_superellipse_p
       << ' ' << config.blending_correction_y_smoothing_passes
       << ' ' << config.state_interpolation_order
       << ' ' << config.derivative_interpolation_order
       << ' ' << config.blending_function_mode
       << ' ' << config.use_sequential_blending
       << ' ' << config.use_c3_repair
       << ' ' << config.c3_apply_x
       << ' ' << config.c3_apply_y
       << ' ' << config.c3_use_interface
       << ' ' << config.c3_buffer_width
       << ' ' << config.c3_boundary_penalty
       << ' ' << config.c3_coefficient_regularization
       << ' ' << config.c3_x_pass_y_smoothing_passes
       << ' ' << config.c3_x_pass_corner_taper
       << ' ' << config.c3_x_alpha_y_smoothing_passes
       << ' ' << config.c3_x_alpha_y_regularization
       << ' ' << config.c3_y_alpha_x_regularization
       << ' ' << config.c3_x_alpha_y_mean_regularization
       << ' ' << config.c3_y_alpha_x_mean_regularization
       << ' ' << config.c3_repair_after_filter
       << ' ' << config.c3_x_transverse_harmonics
       << ' ' << config.c3_y_transverse_harmonics
       << ' ' << config.c3_x_alpha_y_projection_harmonics
       << ' ' << config.c3_y_alpha_x_projection_harmonics
       << ' ' << config.filter_periodic_copies_before_linear_derivatives
       << ' ' << config.linear_derivative_spectral_filter_mode
       << ' ' << config.enable_post_timestep_low_pass_filter
       << ' ' << config.enable_post_rk4_state_filter_and_blend
       << ' ' << config.enable_final_output_blended_derivative_filter
       << ' ' << config.enable_initial_c3_effect_diagnostic
       << ' ' << config.enable_nested_initialization_validation
       << '\n';
}

void readNestedConfig(std::istream& is, crescent::nested::AlgorithmConfig& config)
{
    requireToken(is, "nested_config");
    is >> config.use_crescent_algorithm
       >> config.l1_box_x_lo
       >> config.l1_box_x_hi
       >> config.l1_box_y_lo
       >> config.l1_box_y_hi
       >> config.roi_accurate_fraction
       >> config.roi_spectral_fraction
       >> config.damping_steepness
       >> config.blending_steepness
       >> config.blending_tensor_alpha
       >> config.blending_superellipse_gamma
       >> config.blending_superellipse_p
       >> config.blending_correction_y_smoothing_passes
       >> config.state_interpolation_order
       >> config.derivative_interpolation_order
       >> config.blending_function_mode
       >> config.use_sequential_blending
       >> config.use_c3_repair
       >> config.c3_apply_x
       >> config.c3_apply_y
       >> config.c3_use_interface
       >> config.c3_buffer_width
       >> config.c3_boundary_penalty
       >> config.c3_coefficient_regularization
       >> config.c3_x_pass_y_smoothing_passes
       >> config.c3_x_pass_corner_taper
       >> config.c3_x_alpha_y_smoothing_passes
       >> config.c3_x_alpha_y_regularization
       >> config.c3_y_alpha_x_regularization
       >> config.c3_x_alpha_y_mean_regularization
       >> config.c3_y_alpha_x_mean_regularization
       >> config.c3_repair_after_filter
       >> config.c3_x_transverse_harmonics
       >> config.c3_y_transverse_harmonics
       >> config.c3_x_alpha_y_projection_harmonics
       >> config.c3_y_alpha_x_projection_harmonics
       >> config.filter_periodic_copies_before_linear_derivatives
       >> config.linear_derivative_spectral_filter_mode
       >> config.enable_post_timestep_low_pass_filter
       >> config.enable_post_rk4_state_filter_and_blend
       >> config.enable_final_output_blended_derivative_filter
       >> config.enable_initial_c3_effect_diagnostic
       >> config.enable_nested_initialization_validation;
}

void writeValidStateField(
    const amrex::MultiFab& field,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    const std::string& file_prefix)
{
    amrex::MultiFab valid_state(ba, dm, 1, 0);
    amrex::MultiFab::Copy(valid_state, field, 0, 0, 1, 0);
    amrex::VisMF::Write(valid_state, file_prefix);
}

void readValidStateField(amrex::MultiFab& field, const std::string& file_prefix)
{
    amrex::MultiFab valid_state(field.boxArray(), field.DistributionMap(), 1, 0);
    amrex::VisMF::Read(valid_state, file_prefix);
    amrex::MultiFab::Copy(field, valid_state, 0, 0, 1, 0);
}

std::string checkpointFieldPath(
    int lev,
    const std::string& checkpoint_name,
    const char* field_name)
{
    return amrex::MultiFabFileFullPrefix(
        lev,
        checkpoint_name,
        "Level_",
        field_name);
}

template <typename FieldOperation>
void forEachCheckpointedSurfaceField(
    const crescent::shared::SolutionState& solution,
    int lev,
    FieldOperation&& op)
{
    op("eta", solution.eta[lev]);
    op("phi", solution.phi[lev]);
    op("pressure", solution.pressure[lev]);
}

template <typename FieldOperation>
void forEachCheckpointedSurfaceField(
    crescent::shared::SolutionState& solution,
    int lev,
    FieldOperation&& op)
{
    op("eta", solution.eta[lev]);
    op("phi", solution.phi[lev]);
    op("pressure", solution.pressure[lev]);
}

void writeCheckpointHeader(
    const std::string& checkpoint_name,
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    const auto mesh = reporting_context.meshHierarchy();
    const auto level_context = reporting_context.levelContext(0);
    const auto core = level_context.coreAccess();
    const auto& wave = core.wave;
    const auto& nested_config = reporting_context.nestedConfig();
    const int finest_level = reporting_context.finestLevel();

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string header_file_name = checkpoint_name + "/Header";
    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
    std::ofstream header_file;
    header_file.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    header_file.open(
        header_file_name.c_str(),
        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!header_file.good()) {
        amrex::FileOpenFailed(header_file_name);
    }

    header_file.precision(17);
    header_file << "Checkpoint file for CRESCENT\n";
    header_file << "format_version " << CheckpointFormatVersion << '\n';
    header_file << "finest_level " << finest_level << '\n';
    header_file << "state_vector_size " << runtime.istep.size() << '\n';
    writeVector(header_file, "istep", runtime.istep);
    writeVector(header_file, "nsubsteps", runtime.nsubsteps);
    writeVector(header_file, "t_new", runtime.t_new);
    writeVector(header_file, "t_old", runtime.t_old);
    writeVector(header_file, "dt", runtime.dt);
    header_file << "runtime_scalars"
                << ' ' << runtime.start_time
                << ' ' << runtime.dt_override
                << '\n';
    header_file << "reporting_state"
                << ' ' << reporting.energy_initial_total
                << ' ' << reporting.energy_tracking_initialized
                << ' ' << reporting.omega_a_stage_diagnostics_initialized
                << ' ' << reporting.omega_a_edge_diagnostics_initialized
                << '\n';
    header_file << "wave_spectral"
                << ' ' << wave.spectral.perturbation_order
                << ' ' << wave.spectral.beta1_filter
                << ' ' << wave.spectral.beta2_filter
                << ' ' << wave.spectral.peak_wavenumber_override
                << ' ' << wave.spectral.peak_wavenumber_physical
                << '\n';
    writeNestedConfig(header_file, nested_config);

    for (int lev = 0; lev <= finest_level; ++lev) {
        writeGeometryMetadata(header_file, lev, mesh.geom[lev]);
        header_file << "box_array\n";
        mesh.grids[lev].writeOn(header_file);
        header_file << '\n';
    }
}

void restoreRuntimeState(
    std::istream& header_stream,
    crescent::shared::RuntimeControlState& runtime,
    crescent::shared::WaveRuntimeState& wave,
    crescent::nested::AlgorithmConfig& nested_config,
    int& checkpoint_finest_level,
    amrex::Vector<amrex::BoxArray>& checkpoint_grids,
    const amrex::Vector<amrex::Geometry>& geom,
    int configured_max_level)
{
    std::string title;
    std::getline(header_stream, title);
    if (title != "Checkpoint file for CRESCENT") {
        amrex::Abort("Not a CRESCENT checkpoint: " + title);
    }

    requireToken(header_stream, "format_version");
    int format_version = -1;
    header_stream >> format_version;
    if (format_version != CheckpointFormatVersion) {
        amrex::Abort("Unsupported CRESCENT checkpoint format version");
    }

    requireToken(header_stream, "finest_level");
    header_stream >> checkpoint_finest_level;
    if (checkpoint_finest_level < 0) {
        amrex::Abort("CRESCENT checkpoint finest_level must be non-negative");
    }
    if (checkpoint_finest_level > configured_max_level) {
        amrex::Abort("CRESCENT checkpoint finest_level exceeds configured amr.max_level");
    }
    if (checkpoint_finest_level >= static_cast<int>(geom.size())) {
        amrex::Abort("CRESCENT checkpoint finest_level exceeds available geometry metadata");
    }

    requireToken(header_stream, "state_vector_size");
    int state_vector_size = 0;
    header_stream >> state_vector_size;
    if (state_vector_size <= checkpoint_finest_level) {
        amrex::Abort("CRESCENT checkpoint state vector is smaller than finest_level + 1");
    }
    runtime.resizeLevelState(state_vector_size);
    wave.resizeLevelState(state_vector_size);

    readVector(header_stream, "istep", runtime.istep);
    readVector(header_stream, "nsubsteps", runtime.nsubsteps);
    readVector(header_stream, "t_new", runtime.t_new);
    readVector(header_stream, "t_old", runtime.t_old);
    readVector(header_stream, "dt", runtime.dt);

    requireToken(header_stream, "runtime_scalars");
    header_stream >> runtime.start_time
                  >> runtime.dt_override;

    requireToken(header_stream, "reporting_state");
    header_stream >> runtime.reporting.energy_initial_total
                  >> runtime.reporting.energy_tracking_initialized
                  >> runtime.reporting.omega_a_stage_diagnostics_initialized
                  >> runtime.reporting.omega_a_edge_diagnostics_initialized;

    requireToken(header_stream, "wave_spectral");
    header_stream >> wave.spectral.perturbation_order
                  >> wave.spectral.beta1_filter
                  >> wave.spectral.beta2_filter
                  >> wave.spectral.peak_wavenumber_override
                  >> wave.spectral.peak_wavenumber_physical;

    readNestedConfig(header_stream, nested_config);

    checkpoint_grids.resize(checkpoint_finest_level + 1);
    for (int lev = 0; lev <= checkpoint_finest_level; ++lev) {
        validateGeometryMetadata(header_stream, lev, geom[lev]);
        requireToken(header_stream, "box_array");
        checkpoint_grids[lev].readFrom(header_stream);
    }
}

void readCheckpointHeader(
    const std::string& restart_file,
    crescent::shared::RuntimeControlState& runtime,
    crescent::shared::WaveRuntimeState& wave,
    crescent::nested::AlgorithmConfig& nested_config,
    int& checkpoint_finest_level,
    amrex::Vector<amrex::BoxArray>& checkpoint_grids,
    const amrex::Vector<amrex::Geometry>& geom,
    int configured_max_level)
{
    const std::string header_path = restart_file + "/Header";
    if (!amrex::FileSystem::Exists(header_path)) {
        amrex::Abort("CRESCENT restart checkpoint header does not exist: " + header_path);
    }

    amrex::Vector<char> file_data;
    amrex::ParallelDescriptor::ReadAndBcastFile(header_path, file_data);
    if (file_data.empty()) {
        amrex::Abort("CRESCENT restart checkpoint header is empty: " + header_path);
    }
    std::string header_text(file_data.dataPtr(), file_data.size());
    std::istringstream header_stream(header_text, std::istringstream::in);
    restoreRuntimeState(
        header_stream,
        runtime,
        wave,
        nested_config,
        checkpoint_finest_level,
        checkpoint_grids,
        geom,
        configured_max_level);
}

} // namespace

namespace crescent::solver::checkpoint {

void writeCheckpointFile(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    const auto& runtime = reporting_context.runtime();
    const auto& solution = reporting_context.solutionState();
    const auto mesh = reporting_context.meshHierarchy();
    const int finest_level = reporting_context.finestLevel();
    const std::string checkpoint_name = checkpointName(runtime);

    amrex::Print() << "Writing CRESCENT checkpoint " << checkpoint_name << "\n";

    const int nlevels = finest_level + 1;
    amrex::PreBuildDirectorHierarchy(checkpoint_name, "Level_", nlevels, true);
    writeCheckpointHeader(checkpoint_name, reporting_context);

    for (int lev = 0; lev <= finest_level; ++lev) {
        forEachCheckpointedSurfaceField(
            solution,
            lev,
            [&](const char* field_name, const amrex::MultiFab& field)
            {
                writeValidStateField(
                    field,
                    mesh.grids[lev],
                    mesh.dmap[lev],
                    checkpointFieldPath(lev, checkpoint_name, field_name));
            });
    }
}

void readCheckpointFile(const CheckpointRestartContext& restart_context)
{
    auto solver = restart_context.ownership.solverAccess();
    auto& runtime = solver.core.runtime;
    const std::string restart_file = runtime.checkpoint.restart_file;

    amrex::Print() << "Restarting CRESCENT from checkpoint " << restart_file << "\n";

    int checkpoint_finest_level = -1;
    amrex::Vector<amrex::BoxArray> checkpoint_grids;
    readCheckpointHeader(
        restart_file,
        runtime,
        solver.core.wave,
        solver.nested.config,
        checkpoint_finest_level,
        checkpoint_grids,
        restart_context.geom,
        restart_context.max_level);

    if (crescent::nested::hasActiveNestedL1(
            solver.nested.config,
            restart_context.max_level)) {
        solver.nested.runtime.resize(restart_context.max_level + 1);
    }

    restart_context.define_hierarchy(
        checkpoint_finest_level,
        checkpoint_grids,
        runtime.t_new);

    for (int lev = 0; lev <= checkpoint_finest_level; ++lev) {
        forEachCheckpointedSurfaceField(
            solver.core.solution,
            lev,
            [&](const char* field_name, amrex::MultiFab& field)
            {
                readValidStateField(
                    field,
                    checkpointFieldPath(lev, restart_file, field_name));
            });
    }

    runtime.checkpoint.restart_file = restart_file;
    runtime.checkpoint.restarted = true;

    const auto mesh =
        crescent::shared::makeMeshHierarchyView(
            restart_context.grids,
            restart_context.dmap,
            restart_context.geom,
            restart_context.ref_ratio);
    const auto initialization_context =
        crescent::solver::setup::makeSolverInitializationContext(
            solver.core,
            solver.nested,
            mesh,
            checkpoint_finest_level);
    crescent::solver::checkpoint::rebuildRestartDerivedState(initialization_context);
}

} // namespace crescent::solver::checkpoint

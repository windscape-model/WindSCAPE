#include "CRESCENT_RunLoop.H"

#include <string>

#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include "CRESCENT_Checkpoint.H"
#include "CRESCENT_FreeSurfaceEvolution.H"
#include "CRESCENT_FreeSurfaceEvolutionContext.H"
#include "CRESCENT_RuntimeOutputs.H"
#include "CRESCENT_RuntimeReportingContext.H"
#include "CRESCENT_StepDiagnostics.H"

namespace crescent::solver::run_loop {

namespace {

struct ExecutionState
{
    amrex::Real time;
    int last_plot_step;
    int last_tecplot_step;
    int last_checkpoint_step;
};

std::string runtimeOutputDirectory(const std::string& plot_file)
{
    const std::size_t slash = plot_file.find_last_of("/\\");
    return (slash == std::string::npos) ? std::string() : plot_file.substr(0, slash);
}

ExecutionState makeExecutionState(
    const crescent::solver::FreeSurfaceEvolutionContext& evolution_context)
{
    const auto& runtime = evolution_context.runtime();
    return {
        runtime.t_new[0],
        runtime.checkpoint.restarted ? runtime.istep[0] : 0,
        runtime.checkpoint.restarted ? runtime.istep[0] : 0,
        runtime.istep[0]
    };
}

void writeInitialOutputs(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    if (reporting.plot_int > 0) {
        crescent::solver::runtime_outputs::writePlotFile(reporting_context);
    }
    if (reporting.tecplot_int > 0) {
        crescent::solver::runtime_outputs::writeTecplotFile(reporting_context);
    }
}

void prepareRestartOutputDirectory(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    const auto& runtime = reporting_context.runtime();
    if (!runtime.checkpoint.restarted) {
        return;
    }

    const std::string output_directory =
        runtimeOutputDirectory(runtime.reporting.plot_file);
    if (output_directory.empty()) {
        return;
    }

    if (amrex::ParallelDescriptor::IOProcessor() &&
        !amrex::FileSystem::Exists(output_directory)) {
        if (!amrex::UtilCreateDirectory(output_directory, 0755)) {
            amrex::Abort("Failed to create CRESCENT restart output directory: " + output_directory);
        }
    }
    amrex::ParallelDescriptor::Barrier();
}

void writeScheduledOutputs(
    const crescent::solver::RuntimeReportingContext& reporting_context,
    int step,
    int& last_plot_step,
    int& last_tecplot_step)
{
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    if (reporting.plot_int > 0 &&
        (step + 1) % reporting.plot_int == 0) {
        last_plot_step = step + 1;
        crescent::solver::runtime_outputs::writePlotFile(reporting_context);
    }

    if (reporting.tecplot_int > 0 &&
        (step + 1) % reporting.tecplot_int == 0) {
        last_tecplot_step = step + 1;
        crescent::solver::runtime_outputs::writeTecplotFile(reporting_context);
    }
}

void writeScheduledCheckpoints(
    const crescent::solver::RuntimeReportingContext& reporting_context,
    int step,
    int& last_checkpoint_step)
{
    const auto& checkpoint = reporting_context.runtime().checkpoint;
    if (checkpoint.interval > 0 &&
        (step + 1) % checkpoint.interval == 0) {
        last_checkpoint_step = step + 1;
        crescent::solver::checkpoint::writeCheckpointFile(reporting_context);
    }
}

void writeFinalOutputs(
    const crescent::solver::RuntimeReportingContext& reporting_context,
    int last_plot_step,
    int last_tecplot_step)
{
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    if (reporting.plot_int > 0 &&
        runtime.istep[0] > last_plot_step) {
        crescent::solver::runtime_outputs::writePlotFile(reporting_context);
    }

    if (reporting.tecplot_int > 0 &&
        runtime.istep[0] > last_tecplot_step) {
        crescent::solver::runtime_outputs::writeTecplotFile(reporting_context);
    }
}

void writeFinalCheckpoint(
    const crescent::solver::RuntimeReportingContext& reporting_context,
    int last_checkpoint_step)
{
    const auto& runtime = reporting_context.runtime();
    const auto& checkpoint = runtime.checkpoint;
    if (checkpoint.interval > 0 &&
        runtime.istep[0] > last_checkpoint_step) {
        crescent::solver::checkpoint::writeCheckpointFile(reporting_context);
    }
}

void executeSingleStep(
    crescent::solver::FreeSurfaceEvolutionContext& evolution_context,
    crescent::solver::RuntimeReportingContext& reporting_context,
    int step,
    ExecutionState& state)
{
    const auto& runtime = evolution_context.runtime();

    if (runtime.verbosity >= 2) {
        amrex::Print() << "\n Start timestep " << step + 1
                       << ": time = " << state.time
                       << " (dt = " << runtime.dt[0] << ")" << std::endl;
    }

    crescent::solver::free_surface_evolution::advanceAllLevelsForStep(
        evolution_context,
        state.time);

    if (runtime.verbosity >= 2) {
        amrex::Print() << "   Done timestep " << step + 1 << std::endl;
    }

    crescent::solver::free_surface_evolution::updateTimeStateAfterStep(
        evolution_context,
        state.time);
    crescent::solver::runtime_reporting::recordPostStep(
        reporting_context,
        state.time);
    writeScheduledOutputs(
        reporting_context,
        step,
        state.last_plot_step,
        state.last_tecplot_step);
    writeScheduledCheckpoints(
        reporting_context,
        step,
        state.last_checkpoint_step);
}

} // namespace

void execute(crescent::solver::FreeSurfaceEvolutionContext& evolution_context)
{
    auto& runtime = evolution_context.runtime();
    auto state = makeExecutionState(evolution_context);
    auto reporting_context = evolution_context.reportingContext();

    prepareRestartOutputDirectory(reporting_context);

    if (!runtime.checkpoint.restarted) {
        writeInitialOutputs(reporting_context);
    }
    crescent::solver::runtime_reporting::initialize(
        reporting_context,
        state.time);

    for (int step = runtime.istep[0];
         step < runtime.max_step && state.time < runtime.stop_time;
         ++step) {
        executeSingleStep(
            evolution_context,
            reporting_context,
            step,
            state);

        if (state.time >= runtime.stop_time - 1.0e-6 * runtime.dt[0]) {
            break;
        }
    }

    writeFinalOutputs(
        reporting_context,
        state.last_plot_step,
        state.last_tecplot_step);
    writeFinalCheckpoint(
        reporting_context,
        state.last_checkpoint_step);
}

} // namespace crescent::solver::run_loop

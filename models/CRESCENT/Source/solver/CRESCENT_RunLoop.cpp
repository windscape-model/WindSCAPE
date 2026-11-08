#include "CRESCENT_RunLoop.H"

#include <string>

#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include "CRESCENT_Checkpoint.H"
#include "CRESCENT_FreeSurfaceEvolution.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_RuntimeOutputs.H"
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
    const crescent::shared::SolverContext& context)
{
    const auto& runtime = context.runtime();
    return {
        runtime.t_new[0],
        runtime.checkpoint.restarted ? runtime.istep[0] : 0,
        runtime.checkpoint.restarted ? runtime.istep[0] : 0,
        runtime.istep[0]
    };
}

void writeInitialOutputs(
    const crescent::shared::SolverContext& context)
{
    const auto& runtime = context.runtime();
    const auto& reporting = runtime.reporting;
    if (reporting.plot_int > 0) {
        crescent::solver::runtime_outputs::writePlotFile(context);
    }
    if (reporting.tecplot_int > 0) {
        crescent::solver::runtime_outputs::writeTecplotFile(context);
    }
}

void prepareRestartOutputDirectory(
    const crescent::shared::SolverContext& context)
{
    const auto& runtime = context.runtime();
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
    const crescent::shared::SolverContext& context,
    int step,
    int& last_plot_step,
    int& last_tecplot_step)
{
    const auto& runtime = context.runtime();
    const auto& reporting = runtime.reporting;
    if (reporting.plot_int > 0 &&
        (step + 1) % reporting.plot_int == 0) {
        last_plot_step = step + 1;
        crescent::solver::runtime_outputs::writePlotFile(context);
    }

    if (reporting.tecplot_int > 0 &&
        (step + 1) % reporting.tecplot_int == 0) {
        last_tecplot_step = step + 1;
        crescent::solver::runtime_outputs::writeTecplotFile(context);
    }
}

void writeScheduledCheckpoints(
    const crescent::shared::SolverContext& context,
    int step,
    int& last_checkpoint_step)
{
    const auto& checkpoint = context.runtime().checkpoint;
    if (checkpoint.interval > 0 &&
        (step + 1) % checkpoint.interval == 0) {
        last_checkpoint_step = step + 1;
        crescent::solver::checkpoint::writeCheckpointFile(context);
    }
}

void writeFinalOutputs(
    const crescent::shared::SolverContext& context,
    int last_plot_step,
    int last_tecplot_step)
{
    const auto& runtime = context.runtime();
    const auto& reporting = runtime.reporting;
    if (reporting.plot_int > 0 &&
        runtime.istep[0] > last_plot_step) {
        crescent::solver::runtime_outputs::writePlotFile(context);
    }

    if (reporting.tecplot_int > 0 &&
        runtime.istep[0] > last_tecplot_step) {
        crescent::solver::runtime_outputs::writeTecplotFile(context);
    }
}

void writeFinalCheckpoint(
    const crescent::shared::SolverContext& context,
    int last_checkpoint_step)
{
    const auto& runtime = context.runtime();
    const auto& checkpoint = runtime.checkpoint;
    if (checkpoint.interval > 0 &&
        runtime.istep[0] > last_checkpoint_step) {
        crescent::solver::checkpoint::writeCheckpointFile(context);
    }
}

void executeSingleStep(
    crescent::shared::SolverContext& context,
    int step,
    ExecutionState& state)
{
    const auto& runtime = context.runtime();

    if (runtime.verbosity >= 2) {
        amrex::Print() << "\n Start timestep " << step + 1
                       << ": time = " << state.time
                       << " (dt = " << runtime.dt[0] << ")" << std::endl;
    }

    crescent::solver::free_surface_evolution::advanceAllLevelsForStep(
        context,
        state.time);

    if (runtime.verbosity >= 2) {
        amrex::Print() << "   Done timestep " << step + 1 << std::endl;
    }

    crescent::solver::free_surface_evolution::updateTimeStateAfterStep(
        context,
        state.time);
    crescent::solver::runtime_reporting::recordPostStep(
        context,
        state.time);
    writeScheduledOutputs(
        context,
        step,
        state.last_plot_step,
        state.last_tecplot_step);
    writeScheduledCheckpoints(
        context,
        step,
        state.last_checkpoint_step);
}

} // namespace

void execute(crescent::shared::SolverContext& context)
{
    auto& runtime = context.runtime();
    auto state = makeExecutionState(context);

    prepareRestartOutputDirectory(context);

    // The initial output bundles are populated before the step-0 files.
    crescent::solver::runtime_reporting::initialize(
        context,
        state.time);
    if (!runtime.checkpoint.restarted) {
        writeInitialOutputs(context);
    }

    for (int step = runtime.istep[0];
         step < runtime.max_step && state.time < runtime.stop_time;
         ++step) {
        executeSingleStep(
            context,
            step,
            state);

        if (state.time >= runtime.stop_time - 1.0e-6 * runtime.dt[0]) {
            break;
        }
    }

    writeFinalOutputs(
        context,
        state.last_plot_step,
        state.last_tecplot_step);
    writeFinalCheckpoint(
        context,
        state.last_checkpoint_step);
}

} // namespace crescent::solver::run_loop

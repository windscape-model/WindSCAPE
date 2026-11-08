#include "CRESCENT_StepDiagnostics.H"

#include <cmath>

#include <AMReX_Print.H>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedCoreSnapshots.H"
#include "CRESCENT_NestedEnergy.H"
#include "CRESCENT_NestedTimestep.H"
#include "CRESCENT_PostStep.H"
#include "CRESCENT_SolverContext.H"
#include "CRESCENT_RuntimeDiagnostics.H"

namespace {

struct CoreNestedEnergySnapshot
{
    amrex::Real core_energy;
    amrex::Real nested_energy;
    amrex::Real nested_energy_omega_a;
    amrex::Real core_energy_under_omega_a;
    amrex::Real composite_energy;
};

bool hasActiveNestedDiagnostics(
    const crescent::shared::SolverContext& reporting_context)
{
    return crescent::nested::hasActiveNestedL1(
        reporting_context.nestedConfig(),
        reporting_context.finestLevel());
}

CoreNestedEnergySnapshot computeCoreNestedEnergySnapshot(
    const crescent::shared::SolverContext& reporting_context)
{
    const auto& solution = reporting_context.solutionState();
    // L1 kinetic energy uses the blended derivative bundle's phi_z: it is the
    // vertical velocity consistent with the fields the nested path actually
    // evolves. The standard bundle's phi_z[1] comes from a periodic-FFT DNO
    // applied to the raw non-periodic L1 state and carries boundary Gibbs
    // noise, so it is not a faithful energy measurement.
    const auto& nested_runtime = reporting_context.nestedRuntime();
    const auto& blended_phi_z =
        nested_runtime.coupling_derivatives.l1_blended_derivatives.phi_z[1];
    const auto& nested_derivatives =
        nested_runtime.coupling_derivatives.l1_blended_derivatives;
    const amrex::Box& box_a = nested_runtime.coupling_zones.domains[1].box_A;

    const amrex::Geometry& geom0 = reporting_context.levelGeometry(0);
    const amrex::Geometry& geom1 = reporting_context.levelGeometry(1);

    const amrex::Real core_energy =
        crescent::nested::energy::computeTotalEnergy(
            solution.eta[0],
            solution.phi[0],
            solution.eta_x[0],
            solution.eta_y[0],
            solution.phi_x[0],
            solution.phi_y[0],
            solution.phi_z[0],
            geom0);
    const amrex::Real nested_energy =
        crescent::nested::energy::computeTotalEnergy(
            solution.eta[1],
            solution.phi[1],
            nested_derivatives.eta_x[1],
            nested_derivatives.eta_y[1],
            nested_derivatives.phi_x[1],
            nested_derivatives.phi_y[1],
            blended_phi_z,
            geom1);
    const amrex::Real nested_energy_omega_a =
        crescent::nested::energy::computeTotalEnergyInBox(
            solution.eta[1],
            solution.phi[1],
            nested_derivatives.eta_x[1],
            nested_derivatives.eta_y[1],
            nested_derivatives.phi_x[1],
            nested_derivatives.phi_y[1],
            blended_phi_z,
            geom1,
            box_a);

    // Omega_A footprint in L0 index space (ref ratio from the cell sizes).
    const amrex::IntVect ref_ratio(
        static_cast<int>(std::lround(geom0.CellSize(0) / geom1.CellSize(0))),
        static_cast<int>(std::lround(geom0.CellSize(1) / geom1.CellSize(1))),
        1);
    const amrex::Box box_a_on_core = amrex::coarsen(box_a, ref_ratio);
    const amrex::Real core_energy_under_omega_a =
        crescent::nested::energy::computeTotalEnergyInBox(
            solution.eta[0],
            solution.phi[0],
            solution.eta_x[0],
            solution.eta_y[0],
            solution.phi_x[0],
            solution.phi_y[0],
            solution.phi_z[0],
            geom0,
            box_a_on_core);

    return {
        core_energy,
        nested_energy,
        nested_energy_omega_a,
        core_energy_under_omega_a,
        core_energy - core_energy_under_omega_a + nested_energy_omega_a
    };
}

void writeEnergySnapshot(
    const crescent::shared::SolverContext& reporting_context,
    amrex::Real time)
{
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    if (!hasActiveNestedDiagnostics(reporting_context)) {
        return;
    }

    const CoreNestedEnergySnapshot energy =
        computeCoreNestedEnergySnapshot(reporting_context);
    crescent::nested::energy::writeEnergyToCSV(
        "crescent_energy.csv",
        runtime.istep[0],
        time,
        energy.core_energy,
        energy.nested_energy,
        reporting.energy_tracking_initialized ? reporting.energy_initial_total : 0.0,
        energy.nested_energy_omega_a,
        energy.core_energy_under_omega_a,
        energy.composite_energy);
}

// The blended derivative bundle is only filled by the step loop; at a fresh
// start it still holds its defined-but-unwritten storage. Populate the L0
// snapshot sources and the L1 output bundle once so the step-0 energy row is
// measured the same way as every later row; it runs before the initial
// output. All targets are recomputed from scratch at each step, so this does
// not alter the evolution.
void populateNestedOutputBundlesAtStart(
    const crescent::shared::SolverContext& reporting_context)
{
    // The step-end refresher computes the np1 snapshot fresh from the current
    // L0 state at any step count (the step-start variant only computes at
    // istep 0), and the L1 output recompute consumes only the np1 snapshot.
    auto core_level_context = reporting_context.levelContext(0);
    crescent::nested::snapshots::storeCoreInterpolationSourceAtStepEnd(
        core_level_context.finestLevel(),
        0,
        core_level_context.core(),
        core_level_context.nested(),
        core_level_context.meshHierarchy().geom);

    auto nested_level_context = reporting_context.levelContext(1);
    crescent::nested::free_surface_evolution::recomputeL1OutputDerivatives(
        nested_level_context);

    crescent::solver::poststep::computeVelocityOutputs(core_level_context);
    crescent::solver::poststep::computeVelocityOutputs(nested_level_context);
}

} // namespace

namespace crescent::solver::runtime_reporting {

void initialize(
    crescent::shared::SolverContext& reporting_context,
    amrex::Real time)
{
    auto& runtime = reporting_context.runtime();
    auto& reporting = runtime.reporting;
    if (runtime.checkpoint.restarted) {
        return;
    }
    if (!hasActiveNestedDiagnostics(reporting_context)) {
        return;
    }

    populateNestedOutputBundlesAtStart(reporting_context);

    const CoreNestedEnergySnapshot energy =
        computeCoreNestedEnergySnapshot(reporting_context);
    if (runtime.istep[0] == 0) {
        reporting.energy_initial_total = energy.core_energy + energy.nested_energy;
        reporting.energy_tracking_initialized = true;
    } else {
        // Preserve the current behavior for nonzero-start runs: without an
        // initial energy reference, energy_change_percent stays at zero.
        reporting.energy_initial_total = 0.0;
        reporting.energy_tracking_initialized = false;
    }

    crescent::nested::energy::writeEnergyToCSV(
        "crescent_energy.csv",
        runtime.istep[0],
        time,
        energy.core_energy,
        energy.nested_energy,
        reporting.energy_initial_total,
        energy.nested_energy_omega_a,
        energy.core_energy_under_omega_a,
        energy.composite_energy);
    amrex::Print() << "Energy tracking initialized: E_L0=" << energy.core_energy
                   << ", E_L1=" << energy.nested_energy
                   << ", E_total=" << (energy.core_energy + energy.nested_energy)
                   << std::endl;
}

void recordPostStep(
    crescent::shared::SolverContext& reporting_context,
    amrex::Real time)
{
    const auto& runtime = reporting_context.runtime();
    const auto& reporting = runtime.reporting;
    const bool has_active_nested_diagnostics =
        hasActiveNestedDiagnostics(reporting_context);

    if (has_active_nested_diagnostics) {
        writeEnergySnapshot(reporting_context, time);
    }

    if (reporting.write_omega_a_edge_diagnostics == 1 &&
        has_active_nested_diagnostics) {
        auto level_context = reporting_context.levelContext(1);
        crescent::solver::runtime_diagnostics::writeOmegaAEdgeDiagnostics(
            level_context,
            time);
        crescent::solver::runtime_diagnostics::writeNestedOscillationDiagnostics(
            level_context,
            time);
    }
}

} // namespace crescent::solver::runtime_reporting

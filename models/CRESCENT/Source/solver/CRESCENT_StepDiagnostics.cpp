#include "CRESCENT_StepDiagnostics.H"

#include <AMReX_Print.H>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedEnergy.H"
#include "CRESCENT_RuntimeReportingContext.H"
#include "CRESCENT_RuntimeDiagnostics.H"

namespace {

struct CoreNestedEnergySnapshot
{
    amrex::Real core_energy;
    amrex::Real nested_energy;
};

bool hasActiveNestedDiagnostics(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    return crescent::nested::hasActiveNestedL1(
        reporting_context.nestedConfig(),
        reporting_context.finestLevel());
}

CoreNestedEnergySnapshot computeCoreNestedEnergySnapshot(
    const crescent::solver::RuntimeReportingContext& reporting_context)
{
    const auto& solution = reporting_context.solutionState();
    return {
        crescent::nested::energy::computeTotalEnergy(
            solution.eta[0],
            solution.phi[0],
            solution.phi_z[0],
            reporting_context.levelGeometry(0)),
        crescent::nested::energy::computeTotalEnergy(
            solution.eta[1],
            solution.phi[1],
            solution.phi_z[1],
            reporting_context.levelGeometry(1))
    };
}

void writeEnergySnapshot(
    const crescent::solver::RuntimeReportingContext& reporting_context,
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
        reporting.energy_tracking_initialized ? reporting.energy_initial_total : 0.0);
}

} // namespace

namespace crescent::solver::runtime_reporting {

void initialize(
    crescent::solver::RuntimeReportingContext& reporting_context,
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
        reporting.energy_initial_total);
    amrex::Print() << "Energy tracking initialized: E_L0=" << energy.core_energy
                   << ", E_L1=" << energy.nested_energy
                   << ", E_total=" << (energy.core_energy + energy.nested_energy)
                   << std::endl;
}

void recordPostStep(
    crescent::solver::RuntimeReportingContext& reporting_context,
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
    }
}

} // namespace crescent::solver::runtime_reporting

#ifndef CRESCENT_NESTED_STAGE_DIAGNOSTICS_H
#define CRESCENT_NESTED_STAGE_DIAGNOSTICS_H

#include <AMReX_MultiFab.H>
#include <AMReX_Real.H>

#include "CRESCENT_LevelEvolutionContext.H"
#include "CRESCENT_NestedConfig.H"
#include "solver/CRESCENT_RuntimeDiagnostics.H"

namespace crescent::nested::stage_diagnostics {

using LevelEvolutionContext = crescent::shared::LevelEvolutionContext;
using amrex::MultiFab;
using amrex::Real;

struct RK4StageDiagnosticsContext
{
    int lev;
    int irk;
    Real stage_time;
    Real temporal_weight;
    int verbosity;
    bool write_phi_stage_diagnostic;
    MultiFab* phi_after_c3_snapshot = nullptr;

    MultiFab* phiStageSnapshot()
    {
        return phi_after_c3_snapshot;
    }

    const MultiFab* phiStageSnapshot() const
    {
        return phi_after_c3_snapshot;
    }
};

inline RK4StageDiagnosticsContext makeRK4StageDiagnosticsContext(
    const LevelEvolutionContext& level_context,
    int irk,
    Real stage_time,
    Real temporal_weight)
{
    const int lev = level_context.level();
    const auto& runtime = level_context.runtime();
    const auto& reporting = runtime.reporting;
    const bool write_phi_stage_diagnostic =
        (crescent::nested::isActiveNestedL1Level(
             level_context.nestedConfig(),
             lev) &&
         irk == 0 &&
         reporting.write_omega_a_edge_diagnostics == 1 &&
         runtime.istep[lev] % reporting.omega_a_edge_diagnostic_int == 0);

    RK4StageDiagnosticsContext context {
        lev,
        irk,
        stage_time,
        temporal_weight,
        runtime.verbosity,
        write_phi_stage_diagnostic,
        nullptr
    };

    if (context.write_phi_stage_diagnostic) {
        context.phi_after_c3_snapshot =
            &level_context
                 .nestedRuntime()
                 .runtime_diagnostic_scratch[lev]
                 .phi_stage
                 .phi_stage_snapshot;
    }

    return context;
}

inline void reportPreparedRK4Stage(
    const RK4StageDiagnosticsContext& context,
    const LevelEvolutionContext& level_context,
    const MultiFab& phi_raw,
    const MultiFab& phi_filtered)
{
    if (context.write_phi_stage_diagnostic && context.phiStageSnapshot() != nullptr) {
        crescent::solver::runtime_diagnostics::writeOmegaAEdgePhiStageDiagnostics(
            level_context,
            context.irk,
            context.stage_time,
            phi_raw,
            *context.phiStageSnapshot(),
            phi_filtered);
    }
}

} // namespace crescent::nested::stage_diagnostics

#endif // CRESCENT_NESTED_STAGE_DIAGNOSTICS_H

#ifndef CRESCENT_RUNTIME_DIAGNOSTICS_H
#define CRESCENT_RUNTIME_DIAGNOSTICS_H

#include <AMReX_MultiFab.H>
#include <AMReX_Real.H>

namespace crescent::shared {
class LevelEvolutionContext;
}

namespace crescent::solver::runtime_diagnostics {

void writeOmegaAEdgeDiagnostics(
    crescent::shared::LevelEvolutionContext& level_context,
    amrex::Real time_next);

void writeOmegaAEdgePhiStageDiagnostics(
    const crescent::shared::LevelEvolutionContext& level_context,
    int irk,
    amrex::Real time_stage,
    const amrex::MultiFab& phi_raw,
    const amrex::MultiFab& phi_c3,
    const amrex::MultiFab& phi_filtered);

} // namespace crescent::solver::runtime_diagnostics

#endif // CRESCENT_RUNTIME_DIAGNOSTICS_H

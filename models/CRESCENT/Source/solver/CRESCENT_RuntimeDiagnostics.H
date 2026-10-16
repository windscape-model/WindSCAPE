#ifndef CRESCENT_RUNTIME_DIAGNOSTICS_H
#define CRESCENT_RUNTIME_DIAGNOSTICS_H

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace crescent::shared {
class LevelContext;
}

namespace crescent::solver::runtime_diagnostics {

void writeOmegaAEdgeDiagnostics(
    crescent::shared::LevelContext& level_context,
    amrex::Real time_next);

// Per-step scalar oscillation + periodicity-seam metrics for the nested L1
// level, written to nested_oscillation_diagnostic.csv. Second differences of
// eta detect grid-scale/high-frequency packets in and around Omega_A; the
// seam columns measure how periodic the C3-repaired scratch copies actually
// are (they are the fields every L1 FFT consumes).
void writeNestedOscillationDiagnostics(
    crescent::shared::LevelContext& level_context,
    amrex::Real time_next);

} // namespace crescent::solver::runtime_diagnostics

#endif // CRESCENT_RUNTIME_DIAGNOSTICS_H

#pragma once

#include "CRESCENT_SolverContext.H"

namespace crescent::nested::setup {

// Fresh run: coupling zones and child storage, the child state interpolated
// from the parent, zero child output fields.
void initializeNestedL1Level(
    const crescent::shared::LevelContext& level_context);

// Coupling zones, the coupling and C3 workspaces and the diagnostic scratch
// of the child level (also used on restart).
void initializeDomainsAndStorage(
    const crescent::shared::LevelContext& level_context);

// The live child state is a non-periodic window of the parent: its output
// bundle is assembled later by the coupling, never by an FFT of the state.
void zeroDerivedFields(
    const crescent::shared::LevelContext& level_context);

} // namespace crescent::nested::setup

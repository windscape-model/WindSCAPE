#pragma once

#include "CRESCENT_ConstructionContext.H"
#include "CRESCENT_MeshHierarchyView.H"

namespace crescent::shared {

struct RuntimeControlState;

} // namespace crescent::shared

namespace crescent::solver::setup {
void initializeLevelTimeState(
    const ConstructionContext& construction_context);

void initializeConfiguredConstructionState(
    const ConstructionContext& construction_context);

void printRefinementSummary(
    const shared::RuntimeControlState& runtime,
    shared::MeshHierarchyView mesh,
    int finest_level);

void initializeScratchLevelStorage(
    const ScratchLevelContext& scratch_level_context);

} // namespace crescent::solver::setup

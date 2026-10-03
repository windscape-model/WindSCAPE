#pragma once

namespace crescent::solver::setup {
class SolverInitializationContext;
}

namespace crescent::solver::checkpoint {

void rebuildRestartDerivedState(
    const crescent::solver::setup::SolverInitializationContext& initialization_context);

} // namespace crescent::solver::checkpoint

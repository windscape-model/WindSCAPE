#pragma once

#include <string>

namespace crescent {
struct SolverState;
}

namespace crescent::solver::setup {

// Reads and validates every input of the solver_pp_prefix ParmParse table
// ("crescent" by default; and max_step, start_time and stop_time) into the
// solver state.
void readParameters(
    crescent::SolverState& state,
    int max_level,
    const std::string& solver_pp_prefix = "crescent");

} // namespace crescent::solver::setup

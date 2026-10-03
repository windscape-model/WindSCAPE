#ifndef CRESCENT_CORE_TIMESTEP_LOGGING_H
#define CRESCENT_CORE_TIMESTEP_LOGGING_H

#include <AMReX_Print.H>

#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"

namespace crescent::core::logging {

inline void logRKStageHeader(
    const crescent::shared::RuntimeControlState& runtime_state,
    int rk_stage)
{
    if (runtime_state.verbosity >= 2 && rk_stage > 0) {
        amrex::Print() << "\n--- RK4 stage " << rk_stage << " ---" << std::endl;
    }
}

inline void logPostAdvanceSurfaceState(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev,
    int rk_stage)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "After advanceSurfaceRK4(" << rk_stage << "): eta min/max = "
                   << solution_state.eta[lev].min(0) << " / "
                   << solution_state.eta[lev].max(0) << std::endl;
    if (rk_stage > 0) {
        amrex::Print() << "After advanceSurfaceRK4(" << rk_stage
                       << "): phi min/max = " << solution_state.phi[lev].min(0)
                       << " / " << solution_state.phi[lev].max(0) << std::endl;
    }
}

inline void logPreSurfaceRhsDiagnostics(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev,
    int rk_stage)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    if (rk_stage == 0) {
        amrex::Print() << "Before computeSurfaceRHS(0): eta_x min/max = "
                       << solution_state.eta_x[lev].min(0) << " / "
                       << solution_state.eta_x[lev].max(0) << std::endl;
        amrex::Print() << "Before computeSurfaceRHS(0): phi_z min/max = "
                       << solution_state.phi_z[lev].min(0) << " / "
                       << solution_state.phi_z[lev].max(0) << std::endl;
        return;
    }

    amrex::Print() << "After computeVerticalVelocity(" << rk_stage
                   << "): phi_z min/max = " << solution_state.phi_z[lev].min(0)
                   << " / " << solution_state.phi_z[lev].max(0) << std::endl;
}

inline void logPostSurfaceRhsDiagnostics(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev,
    int rk_stage)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "After computeSurfaceRHS(" << rk_stage
                   << "): rhs_eta[" << rk_stage << "] min/max = "
                   << solution_state.rhs_eta[lev].min(rk_stage) << " / "
                   << solution_state.rhs_eta[lev].max(rk_stage) << std::endl;
    amrex::Print() << "After computeSurfaceRHS(" << rk_stage
                   << "): rhs_phi[" << rk_stage << "] min/max = "
                   << solution_state.rhs_phi[lev].min(rk_stage) << " / "
                   << solution_state.rhs_phi[lev].max(rk_stage) << std::endl;
}

} // namespace crescent::core::logging

#endif // CRESCENT_CORE_TIMESTEP_LOGGING_H

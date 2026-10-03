#ifndef CRESCENT_TIMESTEP_LOGGING_H
#define CRESCENT_TIMESTEP_LOGGING_H

#include <algorithm>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"

namespace crescent::solver::logging {

inline void printMultiFabLeadingRowSample(const amrex::MultiFab& mf, const char* name)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto mf_arr = mf.const_array(mfi);
        const int sample_count = std::min(5, bx.length(0));
        const int i_begin = bx.smallEnd(0);
        const int j = bx.smallEnd(1);
        const int k = bx.smallEnd(2);
        amrex::Print() << "  " << name << "(i=" << i_begin << ":"
                       << (i_begin + sample_count - 1)
                       << ", j=" << j << ", k=" << k << ") = ";
        for (int n = 0; n < sample_count; ++n) {
            amrex::Print() << mf_arr(i_begin + n, j, k);
            if (n < sample_count - 1) {
                amrex::Print() << ", ";
            }
        }
        amrex::Print() << std::endl;
        break;
    }
}

inline void logAdvanceWaveStateSnapshot(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev,
    const char* label)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << label << ": eta min/max = "
                   << solution_state.eta[lev].min(0) << " / "
                   << solution_state.eta[lev].max(0) << std::endl;
    amrex::Print() << label << ": phi min/max = "
                   << solution_state.phi[lev].min(0) << " / "
                   << solution_state.phi[lev].max(0) << std::endl;

    printMultiFabLeadingRowSample(solution_state.eta[lev], "eta");
    printMultiFabLeadingRowSample(solution_state.phi[lev], "phi");
}

inline void logAdvanceWavePreFilterDiagnostics(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev,
    amrex::Real beta1,
    amrex::Real beta2)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "Before  (advance_wave):" << std::endl;
    amrex::Print() << "  eta min/max: " << solution_state.eta[lev].min(0)
                   << " / " << solution_state.eta[lev].max(0) << std::endl;
    amrex::Print() << "  phi min/max: " << solution_state.phi[lev].min(0)
                   << " / " << solution_state.phi[lev].max(0) << std::endl;
    amrex::Print() << "  beta1 = " << beta1 << ", beta2 = " << beta2 << std::endl;
}

inline void logAdvanceWavePostFilterDiagnostics(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "After  (advance_wave):" << std::endl;
    amrex::Print() << "  eta min/max: " << solution_state.eta[lev].min(0)
                   << " / " << solution_state.eta[lev].max(0) << std::endl;
    amrex::Print() << "  phi min/max: " << solution_state.phi[lev].min(0)
                   << " / " << solution_state.phi[lev].max(0) << std::endl;

    printMultiFabLeadingRowSample(solution_state.eta[lev], "eta");
    printMultiFabLeadingRowSample(solution_state.phi[lev], "phi");
}

inline void logAdvanceWaveDerivativeRecomputeDiagnostics(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "After runStandardDerivativeDNOComputation: eta_x min/max = "
                   << solution_state.eta_x[lev].min(0) << " / "
                   << solution_state.eta_x[lev].max(0) << std::endl;
    amrex::Print() << "After runStandardDerivativeDNOComputation: eta min/max = "
                   << solution_state.eta[lev].min(0) << " / "
                   << solution_state.eta[lev].max(0) << std::endl;
}

inline void logAdvanceWaveEndBanner(
    const crescent::shared::RuntimeControlState& runtime_state)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "=== ADVANCE_WAVE END ===" << std::endl << std::endl;
}

inline void logFinalWaveTimestepState(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    int lev)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    amrex::Print() << "\nFINAL STATE at end of timestep "
                   << runtime_state.istep[lev] + 1 << ":" << std::endl;
    amrex::Print() << "  eta min/max = " << solution_state.eta[lev].min(0)
                   << " / " << solution_state.eta[lev].max(0) << std::endl;
    amrex::Print() << "  phi min/max = " << solution_state.phi[lev].min(0)
                   << " / " << solution_state.phi[lev].max(0) << std::endl;
    amrex::Print() << "  phi_z min/max = " << solution_state.phi_z[lev].min(0)
                   << " / " << solution_state.phi_z[lev].max(0) << std::endl;
    amrex::Print() << "  u min/max = " << solution_state.u[lev].min(0)
                   << " / " << solution_state.u[lev].max(0) << std::endl;
    amrex::Print() << "  v min/max = " << solution_state.v[lev].min(0)
                   << " / " << solution_state.v[lev].max(0) << std::endl;

    printMultiFabLeadingRowSample(solution_state.eta[lev], "eta");
    printMultiFabLeadingRowSample(solution_state.phi[lev], "phi");
    printMultiFabLeadingRowSample(solution_state.phi_z[lev], "phi_z");
}

} // namespace crescent::solver::logging

#endif // CRESCENT_TIMESTEP_LOGGING_H

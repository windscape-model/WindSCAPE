#ifndef CRESCENT_NESTED_RUNTIME_LOGGING_H
#define CRESCENT_NESTED_RUNTIME_LOGGING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

namespace crescent::nested::runtime_logging {

inline void logNestedL1StageHeader(int irk,
                                   amrex::Real temporal_weight,
                                   int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "\n--- CRESCENT RK4 Stage " << irk
                   << " (c=" << temporal_weight << ") ---" << std::endl;
}

inline void logNestedStateFilterBlendStart(int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "\n=== Filter and Blend State Variables (after RK4) ==="
                   << std::endl;
}

inline void logNestedStateRangePair(const char* label,
                                    const amrex::MultiFab& eta,
                                    const amrex::MultiFab& phi,
                                    int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "  " << label << ": eta range=["
                   << eta.min(0) << ", " << eta.max(0) << "]" << std::endl;
    amrex::Print() << "  " << label << ": phi range=["
                   << phi.min(0) << ", " << phi.max(0) << "]" << std::endl;
}

inline void logNestedStateFilterBlendComplete(
    const amrex::MultiFab& eta,
    const amrex::MultiFab& phi,
    int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "  After blending: eta range=["
                   << eta.min(0) << ", " << eta.max(0) << "]" << std::endl;
    amrex::Print() << "  After blending: phi range=["
                   << phi.min(0) << ", " << phi.max(0) << "]" << std::endl;
    amrex::Print() << "=== Filter and Blend Complete ===" << std::endl;
}

inline void logNestedFinalStateAfterFilter(
    const amrex::MultiFab& eta,
    amrex::Real eta_amp_before_final,
    amrex::Real eta_amp_after_final,
    int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "Final L1 state (after filter): eta min/max = "
                   << eta.min(0) << " / " << eta.max(0)
                   << ", |eta|_inf=" << eta_amp_after_final
                   << " (total change: "
                   << (eta_amp_after_final - eta_amp_before_final)
                   << ")" << std::endl;
}

inline void logNestedL1TimestepStart(const amrex::MultiFab& eta, int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "\n=== CRESCENT L1 RK4 Timestep ===" << std::endl;
    amrex::Print() << "Initial L1 state: eta min/max = "
                   << eta.min(0) << " / " << eta.max(0) << std::endl;
}

inline void logNestedOutputDerivativeRecomputeStart(int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "\n=== Recomputing L1 derivatives from final state for output ==="
                   << std::endl;
}

inline void logNestedOutputDerivativeRecomputeComplete(int verbosity)
{
    if (verbosity < 2) {
        return;
    }

    amrex::Print() << "  Recomputed L1 derivatives from final state ✓"
                   << std::endl;
}

} // namespace crescent::nested::runtime_logging

#endif // CRESCENT_NESTED_RUNTIME_LOGGING_H

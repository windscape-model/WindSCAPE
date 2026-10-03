#ifndef CRESCENT_LEVEL_INIT_LOGGING_H
#define CRESCENT_LEVEL_INIT_LOGGING_H

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_WaveState.H"

namespace crescent::solver::logging {

inline void logLevelInitializationState(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    const crescent::shared::WaveRuntimeState& wave_state,
    int lev,
    const char* label)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    const auto filter_settings = wave_state.spectralFilterSettings();
    amrex::Print() << label << "  (initialization):" << std::endl;
    amrex::Print() << "  eta min/max: " << solution_state.eta[lev].min(0)
                   << " / " << solution_state.eta[lev].max(0) << std::endl;
    amrex::Print() << "  phi min/max: " << solution_state.phi[lev].min(0)
                   << " / " << solution_state.phi[lev].max(0) << std::endl;
    amrex::Print() << "  beta1_filter = " << filter_settings.beta1_filter
                   << ", beta2_filter = " << filter_settings.beta2_filter
                   << std::endl;
    if (wave_state.spectral.peak_wavenumber_override >= 0.0) {
        amrex::Print() << "  peak_wavenumber_override = "
                       << wave_state.spectral.peak_wavenumber_override
                       << " (user-specified)" << std::endl;
    }
}

inline void logWavenumberPowerSamples(
    const crescent::shared::SolutionState& solution_state,
    const crescent::shared::RuntimeControlState& runtime_state,
    const crescent::shared::WaveRuntimeState& wave_state,
    int lev)
{
    if (runtime_state.verbosity < 2) {
        return;
    }

    for (amrex::MFIter mfi(solution_state.wavenumber_pow[lev]); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto wavenumber_arr = solution_state.wavenumber_pow[lev].const_array(mfi);

        for (int mode = 0; mode <= 2; ++mode) {
            const int i = mode;
            constexpr int j = 0;
            constexpr int k = 0;
            if (bx.contains(amrex::IntVect(i, j, k))) {
                amrex::Print() << "  wavenumber_pow mode " << mode
                               << " at y=0: wavenumber_pow[0]="
                               << wavenumber_arr(i, j, k, 0);
                if (wave_state.spectral.perturbation_order > 1) {
                    amrex::Print() << ", wavenumber_pow[1]="
                                   << wavenumber_arr(i, j, k, 1);
                }
                amrex::Print() << std::endl;
            }
        }

        if (lev == 0 &&
            bx.contains(amrex::IntVect(5, 0, 0))) {
            amrex::Real wvn_min = 1e30;
            amrex::Real wvn_max = -1e30;
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            for (int j = lo.y; j <= hi.y; ++j) {
                const amrex::Real val = wavenumber_arr(5, j, 0, 0);
                wvn_min = std::min(wvn_min, val);
                wvn_max = std::max(wvn_max, val);
            }

            amrex::Print() << "  wavenumber_pow[0] at mode 5 across all y: min="
                           << wvn_min << ", max=" << wvn_max
                           << ", Δ=" << (wvn_max - wvn_min) << std::endl;
        }
    }
}

} // namespace crescent::solver::logging

#endif // CRESCENT_LEVEL_INIT_LOGGING_H

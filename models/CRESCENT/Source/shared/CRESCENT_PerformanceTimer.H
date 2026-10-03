#pragma once

#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <string>
#include <string_view>

/**
 * @file CRESCENT_PerformanceTimer.H
 * @brief Lightweight scoped timing diagnostics for selected solver kernels.
 */

namespace crescent::diagnostics {

class ScopedTimer {
public:
    explicit ScopedTimer(std::string_view name, int verbosity = 1)
        : enabled_(verbosity >= 2)
    {
        if (enabled_) {
            name_.assign(name);
            start_ = amrex::second();
            amrex::Print() << "Starting timer: " << name_ << "\n";
        }
    }

    ~ScopedTimer() {
        if (enabled_) {
            const amrex::Real elapsed_time = elapsed();
            amrex::Print() << "Timer " << name_ << ": "
                           << elapsed_time << " seconds\n";
        }
    }

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

    amrex::Real elapsed() const {
        return enabled_ ? amrex::second() - start_ : amrex::Real(0.0);
    }

private:
    std::string name_;
    bool enabled_{false};
    amrex::Real start_{0.0};
};

} // namespace crescent::diagnostics

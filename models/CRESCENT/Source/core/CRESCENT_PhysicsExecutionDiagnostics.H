#ifndef CRESCENT_PHYSICS_EXECUTION_DIAGNOSTICS_H
#define CRESCENT_PHYSICS_EXECUTION_DIAGNOSTICS_H

#include "CRESCENT_DNO_Internal.H"
#include "CRESCENT_PhiExpansionResearchDiagnostics.H"
#include "CRESCENT_VerticalVelocityResearchDiagnostics.H"

namespace crescent::physics::execution_diagnostics {

inline bool shouldLogResearchDiagnostics(int verbosity)
{
    return crescent::physics::research_diagnostics::shouldLogDetailedDiagnostics(
        verbosity);
}

struct PhiExpansionDiagnosticHooks
{
    int verbosity;

    auto initialSpectrumLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](
                   const amrex::cMultiFab& spectral_phi_component) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logSurfacePotentialSpectrum(
                spectral_phi_component,
                local_verbosity);
        };
    }

    auto transferSpectrumLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](const amrex::cMultiFab& spectral_term,
                                             int expansion_order,
                                             int eta_order) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logFirstExpansionTransferSpectrum(
                spectral_term,
                expansion_order,
                eta_order,
                local_verbosity);
        };
    }

    auto physicalTransferLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](const amrex::MultiFab& physical_term,
                                             const amrex::MultiFab& eta_power_storage,
                                             const amrex::GeometryData& geomdata,
                                             int expansion_order,
                                             int eta_order) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logFirstExpansionYUniformity(
                physical_term,
                eta_power_storage,
                geomdata,
                expansion_order,
                eta_order,
                local_verbosity);
        };
    }

    auto makeExecutionHooks() const
    {
        return crescent::dno::detail::makePhiExpansionExecutionHooks(
            initialSpectrumLogger(),
            transferSpectrumLogger(),
            physicalTransferLogger());
    }
};

inline PhiExpansionDiagnosticHooks makePhiExpansionDiagnosticHooks(int verbosity)
{
    return {verbosity};
}

struct VerticalVelocityDiagnosticHooks
{
    int perturbation_order;
    int verbosity;

    void logInputs(const amrex::MultiFab& phi_expansion,
                   const amrex::MultiFab& eta_pow,
                   const amrex::MultiFab& wavenumber_pow) const
    {
        if (!shouldLogResearchDiagnostics(verbosity)) {
            return;
        }

        crescent::physics::research_diagnostics::logVerticalVelocityInputs(
            phi_expansion,
            eta_pow,
            wavenumber_pow,
            perturbation_order,
            verbosity);
    }

    auto expansionCoefficientLogger() const
    {
        const int local_perturbation_order = perturbation_order;
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_perturbation_order, local_verbosity, should_log](
                   const amrex::MultiFab& expansion_accum) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityExpansionCoefficients(
                expansion_accum,
                local_perturbation_order,
                local_verbosity);
        };
    }

    auto computationStartLogger() const
    {
        const int local_perturbation_order = perturbation_order;
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_perturbation_order, local_verbosity, should_log]() {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityComputationStart(
                local_perturbation_order,
                local_verbosity);
        };
    }

    auto nonlinearTermStartLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](int eta_order, int cumulative_component) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityNonlinearTermStart(
                eta_order,
                cumulative_component,
                local_verbosity);
        };
    }

    auto nonlinearTermStateLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](const amrex::MultiFab& term_state,
                                             const amrex::MultiFab& eta_power_storage,
                                             int eta_order) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityTermState(
                term_state,
                eta_power_storage,
                eta_order,
                local_verbosity);
        };
    }

    auto nonlinearAccumulationLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](const amrex::MultiFab& accumulated_phi_z) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityAccumulation(
                accumulated_phi_z,
                local_verbosity);
        };
    }

    auto nonlinearSummaryLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](const amrex::MultiFab& accumulated_phi_z) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityNonlinearSummary(
                accumulated_phi_z,
                local_verbosity);
        };
    }

    auto finalLinearStartLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log]() {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityFinalLinearTermStart(
                local_verbosity);
        };
    }

    auto finalLinearPreAccumulationLogger() const
    {
        const int local_verbosity = verbosity;
        const bool should_log = shouldLogResearchDiagnostics(local_verbosity);
        return [local_verbosity, should_log](const amrex::MultiFab& term_state) {
            if (!should_log) {
                return;
            }
            crescent::physics::research_diagnostics::logVerticalVelocityFinalLinearContribution(
                term_state,
                local_verbosity);
        };
    }

    auto finalLinearPostAccumulationLogger() const
    {
        return [](const amrex::MultiFab&, const amrex::MultiFab&) {};
    }

    auto makeExecutionHooks() const
    {
        return crescent::dno::detail::makeVerticalVelocityExecutionHooks(
            expansionCoefficientLogger(),
            computationStartLogger(),
            nonlinearTermStartLogger(),
            nonlinearTermStateLogger(),
            nonlinearAccumulationLogger(),
            nonlinearSummaryLogger(),
            finalLinearStartLogger(),
            finalLinearPreAccumulationLogger(),
            finalLinearPostAccumulationLogger());
    }

    void logPreDealiasing(const amrex::MultiFab& phi_z) const
    {
        if (!shouldLogResearchDiagnostics(verbosity)) {
            return;
        }

        crescent::physics::research_diagnostics::logVerticalVelocityPreDealiasing(
            phi_z,
            verbosity);
    }

    void logPostDealiasing(const amrex::MultiFab& phi_z) const
    {
        if (!shouldLogResearchDiagnostics(verbosity)) {
            return;
        }

        crescent::physics::research_diagnostics::logVerticalVelocityPostDealiasing(
            phi_z,
            verbosity);
    }
};

inline VerticalVelocityDiagnosticHooks makeVerticalVelocityDiagnosticHooks(
    int perturbation_order,
    int verbosity)
{
    return {perturbation_order, verbosity};
}

} // namespace crescent::physics::execution_diagnostics

#endif // CRESCENT_PHYSICS_EXECUTION_DIAGNOSTICS_H

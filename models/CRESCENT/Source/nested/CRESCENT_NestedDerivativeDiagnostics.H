#ifndef CRESCENT_NESTED_DERIVATIVE_DIAGNOSTICS_H
#define CRESCENT_NESTED_DERIVATIVE_DIAGNOSTICS_H

#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_DNO_Internal.H"

namespace crescent::nested::diagnostics {

inline bool shouldLogDetailedNestedDerivativeDiagnostics(int verbosity)
{
    return verbosity >= 2;
}

inline bool shouldLogNestedSurfaceDerivativeDiagnostics(int verbosity)
{
    return verbosity >= 2;
}

inline void logSurfaceDerivativeComputationStart(
    int ni,
    int nj,
    amrex::Real Lx,
    amrex::Real Ly,
    amrex::Real dk_x,
    amrex::Real dk_y,
    int spectral_filter_mode,
    const amrex::Box& box_S,
    int verbosity)
{
    if (!shouldLogNestedSurfaceDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\n=== CRESCENT Surface Derivatives (L1 FFT) ===" << std::endl;
    amrex::Print() << "  L1 dimensions: " << ni << " x " << nj << std::endl;
    amrex::Print() << "  L1 physical size: Lx=" << Lx << " m, Ly=" << Ly << " m" << std::endl;
    amrex::Print() << "  L1 wavenumbers: dk_x=" << dk_x
                   << ", dk_y=" << dk_y << " rad/m" << std::endl;
    amrex::Print() << "  Linear derivative spectral filter mode: "
                   << spectral_filter_mode << std::endl;
    amrex::Print() << "  Ω_S box: " << box_S
                   << " (derivatives extracted from this region)" << std::endl;
    amrex::Print() << "  Ω_S box (for extraction): " << box_S << std::endl;
}

inline void logSurfaceDerivativeRanges(
    const char* field_name,
    const amrex::MultiFab& field_x,
    const amrex::MultiFab& field_y,
    int verbosity)
{
    if (!shouldLogNestedSurfaceDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "  " << field_name << "_x range: ["
                   << field_x.min(0) << ", " << field_x.max(0) << "]" << std::endl;
    amrex::Print() << "  " << field_name << "_y range: ["
                   << field_y.min(0) << ", " << field_y.max(0) << "]" << std::endl;
}

inline void logSurfaceDerivativeComputationEnd(int verbosity)
{
    if (!shouldLogNestedSurfaceDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "=== CRESCENT Surface Derivatives L1 complete ===" << std::endl;
}

inline void logNestedVerticalVelocityInputs(
    const amrex::MultiFab& phi_expansion,
    const amrex::MultiFab& wavenumber_pow,
    int perturbation_order,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\n=== computeVerticalVelocity_L1 ===" << std::endl;
    amrex::Print() << "MATHEMATICAL EXPRESSION:" << std::endl;
    amrex::Print() << "  ∂φ/∂z|_{z=η} = Σ_{k=0}^{M-1} η^k/k! * |∇|^{k+1} * φ_k"
                   << std::endl;
    amrex::Print() << "  where φ_k are horizontal expansion terms" << std::endl;
    amrex::Print() << "  Testing with actual Stokes wave η (full nonlinear DNO)"
                   << std::endl;
    amrex::Print() << std::endl;
    crescent::dno::detail::printLeadingComponentRanges(
        "phi_expansion", phi_expansion, perturbation_order);
    crescent::dno::detail::printLeadingComponentRanges(
        "wavenumber_pow", wavenumber_pow, perturbation_order);
}

inline void logNestedVerticalVelocityExpansionCoefficients(
    const amrex::MultiFab& pk,
    int perturbation_order,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nBuilt pk coefficients for perturbation expansion:"
                   << std::endl;
    crescent::dno::detail::printLeadingComponentRanges("pk", pk, perturbation_order);
}

inline void logNestedVerticalVelocityComputationStart(
    int perturbation_order,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nComputing phi_z using " << perturbation_order
                   << "-order expansion:" << std::endl;
}

inline void logNestedVerticalVelocityNonlinearTermStart(
    int kw,
    int k1,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "    kw=" << kw << ", k1=" << k1
                   << ", using pk[" << k1 << "], eta_pow[" << kw << "]"
                   << std::endl;
}

inline void logNestedVerticalVelocityTermState(
    const amrex::MultiFab& tmp,
    const amrex::MultiFab& eta_pow,
    int kw,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "      tmp min/max after forwardThenBackward = "
                   << tmp.min(0) << " / " << tmp.max(0) << std::endl;
    amrex::Print() << "      eta_pow[" << kw << "] max = "
                   << eta_pow.max(kw) << std::endl;
}

inline void logNestedVerticalVelocityNonlinearSummary(
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nAfter nonlinear terms (kw loop):" << std::endl;
    amrex::Print() << "  phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
}

inline void logNestedVerticalVelocityFinalLinearTermStart(int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nComputing final linear term (wvn_power=1)..."
                   << std::endl;
}

inline void logNestedVerticalVelocityFinalContribution(
    const amrex::MultiFab& tmp,
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "  Final term contribution tmp min/max = "
                   << tmp.min(0) << " / " << tmp.max(0) << std::endl;
    amrex::Print() << "  After final linear term: phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
}

inline void logNestedVerticalVelocitySuccess(
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedNestedDerivativeDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "phi_z computed successfully (range: ["
                   << phi_z.min(0) << ", " << phi_z.max(0) << "])"
                   << std::endl;
}

struct VerticalVelocityDiagnosticHooks
{
    int perturbation_order;
    int verbosity;

    auto makeExecutionHooks() const
    {
        const int local_perturbation_order = perturbation_order;
        const int local_verbosity = verbosity;

        return crescent::dno::detail::makeVerticalVelocityExecutionHooks(
            [local_perturbation_order, local_verbosity](
                const amrex::MultiFab& expansion_accum) {
                logNestedVerticalVelocityExpansionCoefficients(
                    expansion_accum,
                    local_perturbation_order,
                    local_verbosity);
            },
            [local_perturbation_order, local_verbosity]() {
                logNestedVerticalVelocityComputationStart(
                    local_perturbation_order,
                    local_verbosity);
            },
            [local_verbosity](int eta_order, int cumulative_component) {
                logNestedVerticalVelocityNonlinearTermStart(
                    eta_order,
                    cumulative_component,
                    local_verbosity);
            },
            [local_verbosity](const amrex::MultiFab& term_state,
                              const amrex::MultiFab& eta_power_storage,
                              int eta_order) {
                logNestedVerticalVelocityTermState(
                    term_state,
                    eta_power_storage,
                    eta_order,
                    local_verbosity);
            },
            [](const amrex::MultiFab&) {},
            [local_verbosity](const amrex::MultiFab& accumulated_phi_z) {
                logNestedVerticalVelocityNonlinearSummary(
                    accumulated_phi_z,
                    local_verbosity);
            },
            [local_verbosity]() {
                logNestedVerticalVelocityFinalLinearTermStart(
                    local_verbosity);
            },
            [](const amrex::MultiFab&) {},
            [local_verbosity](const amrex::MultiFab& term_state,
                              const amrex::MultiFab& accumulated_phi_z) {
                logNestedVerticalVelocityFinalContribution(
                    term_state,
                    accumulated_phi_z,
                    local_verbosity);
            });
    }

    void logInputs(const amrex::MultiFab& phi_expansion,
                   const amrex::MultiFab& wavenumber_pow) const
    {
        logNestedVerticalVelocityInputs(
            phi_expansion,
            wavenumber_pow,
            perturbation_order,
            verbosity);
    }

    void logSuccess(const amrex::MultiFab& phi_z) const
    {
        logNestedVerticalVelocitySuccess(phi_z, verbosity);
    }
};

inline VerticalVelocityDiagnosticHooks makeVerticalVelocityDiagnosticHooks(
    int perturbation_order,
    int verbosity)
{
    return {perturbation_order, verbosity};
}

} // namespace crescent::nested::diagnostics

#endif // CRESCENT_NESTED_DERIVATIVE_DIAGNOSTICS_H

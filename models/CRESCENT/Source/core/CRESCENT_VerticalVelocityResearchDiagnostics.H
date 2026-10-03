#ifndef CRESCENT_VERTICAL_VELOCITY_RESEARCH_DIAGNOSTICS_H
#define CRESCENT_VERTICAL_VELOCITY_RESEARCH_DIAGNOSTICS_H

#include "CRESCENT_DNO_Internal.H"
#include "CRESCENT_PhysicsResearchDiagnostics.H"

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

namespace crescent::physics::research_diagnostics {

inline void logVerticalVelocityInputs(
    const amrex::MultiFab& phi_expansion,
    const amrex::MultiFab& eta_pow,
    const amrex::MultiFab& wavenumber_pow,
    int perturbation_order,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\n=== WSURF INPUTS ===" << std::endl;
    crescent::dno::detail::printLeadingComponentRanges(
        "phi_expansion", phi_expansion, perturbation_order);
    crescent::dno::detail::printLeadingComponentRanges(
        "eta_pow", eta_pow, perturbation_order);
    crescent::dno::detail::printLeadingComponentRanges(
        "wavenumber_pow", wavenumber_pow, perturbation_order);

    for (amrex::MFIter mfi(phi_expansion); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto phi_arr = phi_expansion.const_array(mfi);
        auto eta_pow_arr = eta_pow.const_array(mfi);
        auto wavenumber_arr = wavenumber_pow.const_array(mfi);
        const int sample_count = std::min(5, bx.length(0));
        const int i_begin = bx.smallEnd(0);
        const int j = bx.smallEnd(1);
        const int k = bx.smallEnd(2);
        amrex::Print() << "First " << sample_count
                       << " elements on valid-box low-y row:" << std::endl;
        printLeadingComponentValuesOnRow(
            "phi_expansion", phi_arr, i_begin, j, k, sample_count, perturbation_order);
        printLeadingComponentValuesOnRow(
            "eta_pow", eta_pow_arr, i_begin, j, k, sample_count, perturbation_order);
        printLeadingComponentValuesOnRow(
            "wavenumber_pow", wavenumber_arr, i_begin, j, k, sample_count, perturbation_order);
        break;
    }
}

inline void logVerticalVelocityExpansionCoefficients(
    const amrex::MultiFab& pk,
    int perturbation_order,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nBuilt pk coefficients for perturbation expansion:"
                   << std::endl;
    crescent::dno::detail::printLeadingComponentRanges("pk", pk, perturbation_order);

    for (amrex::MFIter mfi(pk); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto pk_arr = pk.const_array(mfi);
        const int sample_count = std::min(5, bx.length(0));
        printLeadingComponentValuesOnRow(
            "pk",
            pk_arr,
            bx.smallEnd(0),
            bx.smallEnd(1),
            bx.smallEnd(2),
            sample_count,
            perturbation_order);
        break;
    }
}

inline void logVerticalVelocityComputationStart(
    int perturbation_order,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nComputing phi_z using " << perturbation_order
                   << "-order expansion:" << std::endl;
}

inline void logVerticalVelocityNonlinearTermStart(
    int kw,
    int k1,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "    kw=" << kw << ", k1=" << k1
                   << ", using pk[" << k1 << "], wavenumber_pow[" << (kw + 1)
                   << "], eta_pow[" << kw << "]" << std::endl;
}

inline void logVerticalVelocityTermState(
    const amrex::MultiFab& tmp,
    const amrex::MultiFab& eta_pow,
    int kw,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "      tmp min/max after forwardThenBackward = "
                   << tmp.min(0) << " / " << tmp.max(0) << std::endl;
    amrex::Print() << "      eta_pow[" << kw << "] max = "
                   << eta_pow.max(kw) << std::endl;
}

inline void logVerticalVelocityAccumulation(
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "      After accumulation: phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
}

inline void logVerticalVelocityNonlinearSummary(
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nAfter nonlinear terms (kw loop):" << std::endl;
    amrex::Print() << "  phi_z min/max = " << phi_z.min(0)
                   << " / " << phi_z.max(0) << std::endl;
}

inline void logVerticalVelocityFinalLinearTermStart(int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "\nComputing final linear term (wvn_power=1)..."
                   << std::endl;
}

inline void logVerticalVelocityFinalLinearContribution(
    const amrex::MultiFab& tmp,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "  Final term contribution tmp min/max = "
                   << tmp.min(0) << " / " << tmp.max(0) << std::endl;
}

inline void logVerticalVelocityPreDealiasing(
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "  After adding final term: phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
    amrex::Print() << "  Before dealiasing: phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
}

inline void logVerticalVelocityPostDealiasing(
    const amrex::MultiFab& phi_z,
    int verbosity)
{
    if (!shouldLogDetailedDiagnostics(verbosity)) {
        return;
    }

    amrex::Print() << "  After dealiasing: phi_z min/max = "
                   << phi_z.min(0) << " / " << phi_z.max(0) << std::endl;
}

} // namespace crescent::physics::research_diagnostics

#endif // CRESCENT_VERTICAL_VELOCITY_RESEARCH_DIAGNOSTICS_H

#pragma once

#include <AMReX_MultiFab.H>

namespace crescent::dno {

/**
 * @brief Non-owning workspace for Dirichlet-to-Neumann Operator (DNO)
 *
 * This bundles the scratch arrays already allocated by the solver and keeps
 * the production DNO call boundary compact without taking memory ownership.
 */
struct WorkspaceView
{
    amrex::MultiFab& eta_pow;
    amrex::MultiFab& phi_expansion;

    amrex::MultiFab& expansion_deriv;
    amrex::cMultiFab& deriv_spectral;

    amrex::MultiFab& expansion_accum;
    amrex::MultiFab& spectral_work;

    amrex::cMultiFab& phi_expansion_spectrum;

    static WorkspaceView fromExisting(
        amrex::MultiFab& eta_pow_,
        amrex::MultiFab& phi_expansion_,
        amrex::MultiFab& expansion_deriv_,
        amrex::cMultiFab& deriv_spectral_,
        amrex::MultiFab& expansion_accum_,
        amrex::MultiFab& spectral_work_,
        amrex::cMultiFab& phi_expansion_spectrum_)
    {
        return WorkspaceView(
            eta_pow_,
            phi_expansion_,
            expansion_deriv_,
            deriv_spectral_,
            expansion_accum_,
            spectral_work_,
            phi_expansion_spectrum_);
    }

private:
    WorkspaceView(
        amrex::MultiFab& eta_pow_,
        amrex::MultiFab& phi_expansion_,
        amrex::MultiFab& expansion_deriv_,
        amrex::cMultiFab& deriv_spectral_,
        amrex::MultiFab& expansion_accum_,
        amrex::MultiFab& spectral_work_,
        amrex::cMultiFab& phi_expansion_spectrum_)
        : eta_pow(eta_pow_),
          phi_expansion(phi_expansion_),
          expansion_deriv(expansion_deriv_),
          deriv_spectral(deriv_spectral_),
          expansion_accum(expansion_accum_),
          spectral_work(spectral_work_),
          phi_expansion_spectrum(phi_expansion_spectrum_)
    {}
};

struct PhiExpansionWorkspaceView
{
    amrex::MultiFab& phi_expansion;
    amrex::MultiFab& eta_pow;
    amrex::MultiFab& expansion_deriv;
    amrex::cMultiFab& deriv_spectral;
    amrex::cMultiFab& phi_expansion_spectrum;
};

inline PhiExpansionWorkspaceView makePhiExpansionWorkspaceView(
    WorkspaceView& workspace)
{
    return {
        workspace.phi_expansion,
        workspace.eta_pow,
        workspace.expansion_deriv,
        workspace.deriv_spectral,
        workspace.phi_expansion_spectrum
    };
}

struct VerticalVelocityWorkspaceView
{
    const amrex::MultiFab& phi_expansion;
    const amrex::MultiFab& eta_pow;
    amrex::MultiFab& expansion_accum;
    amrex::MultiFab& spectral_work;
};

inline VerticalVelocityWorkspaceView makeVerticalVelocityWorkspaceView(
    WorkspaceView& workspace)
{
    return {
        workspace.phi_expansion,
        workspace.eta_pow,
        workspace.expansion_accum,
        workspace.spectral_work
    };
}

} // namespace crescent::dno

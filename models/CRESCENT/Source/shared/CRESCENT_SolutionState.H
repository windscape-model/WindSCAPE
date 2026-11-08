#ifndef CRESCENT_SOLUTION_STATE_H
#define CRESCENT_SOLUTION_STATE_H

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "CRESCENT_DerivativeViews.H"

namespace crescent::shared {

struct SolutionState
{
    amrex::Vector<amrex::MultiFab> eta;
    amrex::Vector<amrex::MultiFab> eta_x;
    amrex::Vector<amrex::MultiFab> eta_y;
    amrex::Vector<amrex::MultiFab> phi;
    amrex::Vector<amrex::MultiFab> phi_x;
    amrex::Vector<amrex::MultiFab> phi_y;
    amrex::Vector<amrex::MultiFab> rhs_eta;
    amrex::Vector<amrex::MultiFab> rhs_phi;
    amrex::Vector<amrex::MultiFab> u;
    amrex::Vector<amrex::MultiFab> v;
    amrex::Vector<amrex::MultiFab> phi_z;
    amrex::Vector<amrex::MultiFab> pressure;
    amrex::Vector<amrex::MultiFab> wavenumber_pow;
    amrex::Vector<amrex::MultiFab> phi_expansion;
    amrex::Vector<amrex::MultiFab> eta_pow;
    amrex::Vector<amrex::MultiFab> eta_saved;
    amrex::Vector<amrex::MultiFab> phi_saved;

    void resize(int nlevs_max)
    {
        eta.resize(nlevs_max);
        eta_x.resize(nlevs_max);
        eta_y.resize(nlevs_max);
        phi.resize(nlevs_max);
        phi_x.resize(nlevs_max);
        phi_y.resize(nlevs_max);
        rhs_eta.resize(nlevs_max);
        rhs_phi.resize(nlevs_max);
        u.resize(nlevs_max);
        v.resize(nlevs_max);
        phi_z.resize(nlevs_max);
        pressure.resize(nlevs_max);
        wavenumber_pow.resize(nlevs_max);
        phi_expansion.resize(nlevs_max);
        eta_pow.resize(nlevs_max);
        eta_saved.resize(nlevs_max);
        phi_saved.resize(nlevs_max);
    }

    void saveSurfaceState(int lev)
    {
        amrex::MultiFab::Copy(eta_saved[lev], eta[lev], 0, 0, 1, 0);
        amrex::MultiFab::Copy(phi_saved[lev], phi[lev], 0, 0, 1, 0);
    }

    void defineScratchLevelStorage(
        int lev,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        int perturbation_order)
    {
        const int scalar_ncomp = 1;
        const int scalar_nghost = (lev == 0) ? 1 : 2;

        eta[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        eta_x[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        eta_y[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        phi[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        phi_x[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        phi_y[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        rhs_eta[lev].define(ba, dm, 4, scalar_nghost);
        rhs_phi[lev].define(ba, dm, 4, scalar_nghost);
        phi_saved[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        eta_saved[lev].define(ba, dm, scalar_ncomp, scalar_nghost);
        pressure[lev].define(ba, dm, scalar_ncomp, scalar_nghost);

        phi_expansion[lev].define(ba, dm, perturbation_order, scalar_nghost);
        eta_pow[lev].define(ba, dm, perturbation_order, scalar_nghost);
        // wavenumber_pow is allocated later on the spectral layout, after FFT setup.

        u[lev].define(ba, dm, 1, 1);
        v[lev].define(ba, dm, 1, 1);
        phi_z[lev].define(ba, dm, 1, 1);

        eta[lev].setVal(0.0);
        eta_x[lev].setVal(0.0);
        eta_y[lev].setVal(0.0);
        phi[lev].setVal(0.0);
        phi_x[lev].setVal(0.0);
        phi_y[lev].setVal(0.0);
        rhs_eta[lev].setVal(0.0);
        rhs_phi[lev].setVal(0.0);
        u[lev].setVal(0.0);
        v[lev].setVal(0.0);
        phi_z[lev].setVal(0.0);
        pressure[lev].setVal(0.0);
        phi_expansion[lev].setVal(0.0);
        eta_pow[lev].setVal(0.0);
        eta_saved[lev].setVal(0.0);
        phi_saved[lev].setVal(0.0);
    }

    DerivativeBundleView derivativeBundleAt(int lev)
    {
        return {eta_x[lev], eta_y[lev], phi_x[lev], phi_y[lev], phi_z[lev]};
    }
};

} // namespace crescent::shared

#endif // CRESCENT_SOLUTION_STATE_H

#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "AMReX_GpuQualifiers.H"

#include "ERF_ProbCommon.H"

/**
 * @brief ERF-CRESCENT Coupling Test Problem
 *
 * Test case for coupled atmosphere-ocean wave interaction. Starts flat (Ampl = 0);
 * the wave forcing is supplied dynamically by CRESCENT through the Coupler (surface
 * pressure out, wave-driven moving terrain in). Migrated to the current ProbBase API.
 */
struct ProbParm : ProbParmDefaults {
    amrex::Real T_0 = 300.0;        // Surface temperature (K)
    amrex::Real Ampl = 0.0;         // Initial wave amplitude (m) - zero for flat start
    amrex::Real wavelength = 100.0; // Wavelength (m) - for compatibility
}; // ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitDensityHSE.H"

    void init_custom_pert(
        const amrex::Box& bx,
        amrex::Array4<amrex::Real const> const& state,
        amrex::Array4<amrex::Real      > const& state_pert,
        amrex::Array4<amrex::Real      > const& r_hse,
        amrex::Array4<amrex::Real      > const& p_hse,
        amrex::Array4<amrex::Real const> const& z_nd,
        amrex::Array4<amrex::Real const> const& z_cc,
        amrex::GeometryData const& geomdata,
        amrex::Array4<amrex::Real const> const& mf_m,
        const SolverChoice& sc,
        const int lev) override;

    void erf_init_rayleigh(
        amrex::Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
        amrex::Geometry const& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_nd,
        amrex::Real zdamp) override;

protected:

private:
    ProbParm parms;
};

#endif

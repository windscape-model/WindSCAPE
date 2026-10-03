#include "ERF_Prob.H"
#include "ERF_EOS.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("Ampl", parms.Ampl);
  pp.query("wavelength", parms.wavelength);
  pp.query("T_0" , parms.T_0);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    const SolverChoice& sc,
    const int /*lev*/)
{
  // Coupled run starts flat (Ampl = 0); the wave forcing comes from CRESCENT via the Coupler.
  // The analytic linear-wave perturbation below is retained for standalone testing and is
  // identically zero when Ampl = 0. Velocity perturbations use the ProbBase default (zero),
  // consistent with the flat start. (Migrated to the current init_custom_pert signature.)
  const bool use_moisture = (sc.moisture_type != MoistureType::None);

  Real H           = geomdata.ProbHi()[2];
  Real Ampl        = parms.Ampl;
  Real wavelength  = parms.wavelength;
  Real kp          = 2.0 * PI / wavelength;
  Real g           = CONST_GRAV;
  Real omega       = std::sqrt(g * kp);

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const auto prob_lo  = geomdata.ProbLo();
      const auto dx       = geomdata.CellSize();
      const Real x = prob_lo[0] + (i + 0.5) * dx[0];

      Real z = z_cc(i,j,k);
      Real z_base = Ampl * std::sin(kp * x);
      z -= z_base;

      Real fac     = std::cosh( kp * (z - H) ) / std::sinh(kp * H);
      Real p_prime = -(Ampl * omega * omega / kp) * fac * std::sin(kp * x) * r_hse(i,j,k);
      Real p_total = p_hse(i,j,k) + p_prime;

      // Define (rho theta) given pprime
      state_pert(i, j, k, RhoTheta_comp) = getRhoThetagivenP(p_total) - getRhoThetagivenP(p_hse(i,j,k));

      // Set scalar = 0 everywhere
      state_pert(i, j, k, RhoScalar_comp) = state_pert(i,j,k,Rho_comp);

      if (use_moisture) {
          state_pert(i, j, k, RhoQ1_comp) = 0.0;
          state_pert(i, j, k, RhoQ2_comp) = 0.0;
      }
  });

  amrex::Gpu::streamSynchronize();
}

void
Problem::erf_init_rayleigh (Vector<amrex::Vector<amrex::Real> >& rayleigh_ptrs,
                            Geometry const& geom,
                            std::unique_ptr<MultiFab>& /*z_phys_nd*/,
                            amrex::Real /*zdamp*/)
{
    const int khi = geom.Domain().bigEnd()[2];
    for (int k = 0; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::ubar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::vbar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::wbar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = 300.0;
    }

    // Damping above k = 60
    for (int k = 60; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::ubar][k]     =   2.0;
        rayleigh_ptrs[Rayleigh::vbar][k]     =   1.0;
        rayleigh_ptrs[Rayleigh::wbar][k]     =   0.0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = 300.0;
   }
}

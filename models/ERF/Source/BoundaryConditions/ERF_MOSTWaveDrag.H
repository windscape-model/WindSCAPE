#ifndef ERF_MOST_WAVE_DRAG_H_
#define ERF_MOST_WAVE_DRAG_H_

#include <cmath>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"
#include "ERF_MovingWave.H"

/**
 * Wave drag: the logarithmic law applied at each point to the wind relative to the orbital
 * velocity u_w = a R omega sin(kx - omega t), w_w = -a R omega cos(kx - omega t) of the moving
 * surface h = a R(t) sin(kx - omega t) (Zhu et al. 2023, Mon. Wea. Rev. 151, Eqs. 26-29):
 *
 *   u_r = u - u_w,  v_r = v,  w_r = w - w_w                   relative velocity
 *   u_par = (u_r + w_r dh/dx) / sqrt((dh/dx)^2 + 1)           tangent to the surface (v_par alike)
 *   tau13 = -rho Cd |u_r| u_par,  tau23 = -rho Cd |u_r| v_par,  |u_r| the speed of (u_r, v_r, w_r)
 *   Cd = [kappa / ln(z_ref / z0)]^2
 *
 * With the stability option, Cd = [kappa / (ln(z_ref / z0) - psi_m(z_ref / L))]^2 and the heat
 * flux -rho u_* theta_* follow from the Monin-Obukhov relations at each point along the wave
 * (Sullivan et al. 2014, J. Atmos. Sci. 71, appendix), with the stability functions of
 * Large et al. (1994). The stresses are at the faces of u and v, with the sign of compute_u_flux.
 * During the ramp, w_w omits the part a R'(t) sin(kx - omega t) of dh/dt.
 */
namespace WaveDrag {

    struct StabParams {
        int         mode;     // 0: kinematic surface heat flux Q, 1: surface temperature theta_s
        amrex::Real zref;     // height of the first cell center above the surface
        amrex::Real z0;       // momentum roughness
        amrex::Real z0h;      // thermal roughness
        amrex::Real Q;        // kinematic surface heat flux [K m/s]
        amrex::Real theta_s;  // surface potential temperature [K]
    };

    // Stability functions of Large et al. (1994)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_m_stab (const amrex::Real zeta)
    {
        if (zeta >= amrex::Real(0.0)) { return -amrex::Real(5.0) * zeta; }
        const amrex::Real x = std::pow(amrex::Real(1.0) - amrex::Real(16.0)*zeta, amrex::Real(0.25));
        return amrex::Real(2.0)*std::log(amrex::Real(0.5)*(amrex::Real(1.0)+x))
             + std::log(amrex::Real(0.5)*(amrex::Real(1.0)+x*x))
             - amrex::Real(2.0)*std::atan(x) + amrex::Real(0.5)*amrex::Real(PI);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real psi_h_stab (const amrex::Real zeta)
    {
        if (zeta >= amrex::Real(0.0)) { return -amrex::Real(5.0) * zeta; }
        const amrex::Real x = std::pow(amrex::Real(1.0) - amrex::Real(16.0)*zeta, amrex::Real(0.25));
        return amrex::Real(2.0)*std::log(amrex::Real(0.5)*(amrex::Real(1.0)+x*x));
    }

    // Fixed-point iteration for u_* and theta_* from the relative speed V0 and the first-cell
    // theta, with zeta = z_ref/L limited to [-10, 2]. den_m = ln(z_ref/z0) - psi_m (at least 1).
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void iterate_ustar_tstar (const amrex::Real& V0, const amrex::Real& theta1,
                              const StabParams& sp,
                              amrex::Real& ustar, amrex::Real& tstar, amrex::Real& den_m)
    {
        const amrex::Real ln_z0  = std::log(sp.zref / sp.z0);
        const amrex::Real ln_z0h = std::log(sp.zref / sp.z0h);
        amrex::Real zeta = amrex::Real(0.0);
        den_m = ln_z0;
        ustar = KAPPA * V0 / ln_z0;
        tstar = amrex::Real(0.0);
        for (int it = 0; it < 8; ++it) {
            const amrex::Real psim  = psi_m_stab(zeta);
            const amrex::Real psih  = psi_h_stab(zeta);
            den_m = amrex::max(ln_z0  - psim, amrex::Real(1.0));
            const amrex::Real den_h = amrex::max(ln_z0h - psih, amrex::Real(1.0));
            ustar = KAPPA * V0 / den_m;
            tstar = (sp.mode == 0) ? -sp.Q / amrex::max(ustar, amrex::Real(1.0e-6))
                                   : KAPPA * (theta1 - sp.theta_s) / den_h;
            if (std::abs(tstar) < amrex::Real(1.0e-10)) { zeta = amrex::Real(0.0); continue; }
            const amrex::Real L = ustar*ustar*theta1 / (KAPPA * CONST_GRAV * tstar);
            zeta = sp.zref / L;
            zeta = amrex::min(amrex::max(zeta, amrex::Real(-10.0)), amrex::Real(2.0));
        }
    }

    // Drag coefficient with the stability correction
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real stab_cd (const amrex::Real& V0, const amrex::Real& theta1, const StabParams& sp)
    {
        amrex::Real ustar, tstar, den_m;
        iterate_ustar_tstar(V0, theta1, sp, ustar, tstar, den_m);
        return (KAPPA/den_m) * (KAPPA/den_m);
    }

    // tau13 at the x-face (i,j,k) at x = xu; with stab, Cd comes from stab_cd
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_drag_u (const int& i, const int& j, const int& k,
                         const amrex::Array4<const amrex::Real>& cons_arr,
                         const amrex::Array4<const amrex::Real>& velx_arr,
                         const amrex::Array4<const amrex::Real>& vely_arr,
                         const amrex::Array4<const amrex::Real>& velz_arr,
                         const amrex::Real& xu, const amrex::Real& time,
                         const MovingWave& wave, amrex::Real Cd,
                         const bool stab, const StabParams& sp)
    {
        const amrex::Real r    = wave.ramp(time);
        const amrex::Real ph   = wave.k * xu - wave.omega * time;
        const amrex::Real u_w  =  wave.a * r * wave.omega * std::sin(ph);
        const amrex::Real w_w  = -wave.a * r * wave.omega * std::cos(ph);
        const amrex::Real dhdx =  wave.a * r * wave.k     * std::cos(ph);

        const amrex::Real uloc = velx_arr(i,j,k);
        const amrex::Real vloc = amrex::Real(0.25) * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                                     + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        const amrex::Real wloc = amrex::Real(0.25) * ( velz_arr(i  ,j,k) + velz_arr(i  ,j,k+1)
                                                     + velz_arr(i-1,j,k) + velz_arr(i-1,j,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;
        const amrex::Real w_r = wloc - w_w;

        const amrex::Real u_par = (u_r + w_r * dhdx) / std::sqrt(dhdx*dhdx + amrex::Real(1.0));
        const amrex::Real V0    = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);
        const amrex::Real rho   = amrex::Real(0.5) * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        if (stab) {
            const amrex::Real theta1 = amrex::Real(0.5) * ( cons_arr(i-1,j,k,RhoTheta_comp)/cons_arr(i-1,j,k,Rho_comp)
                                                          + cons_arr(i  ,j,k,RhoTheta_comp)/cons_arr(i  ,j,k,Rho_comp) );
            Cd = stab_cd(V0, theta1, sp);
        }
        return -rho * Cd * V0 * u_par;
    }

    // tau23 at the y-face (i,j,k) at x = xv; the wave travels in x, so dh/dy = 0
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_drag_v (const int& i, const int& j, const int& k,
                         const amrex::Array4<const amrex::Real>& cons_arr,
                         const amrex::Array4<const amrex::Real>& velx_arr,
                         const amrex::Array4<const amrex::Real>& vely_arr,
                         const amrex::Array4<const amrex::Real>& velz_arr,
                         const amrex::Real& xv, const amrex::Real& time,
                         const MovingWave& wave, amrex::Real Cd,
                         const bool stab, const StabParams& sp)
    {
        const amrex::Real r    = wave.ramp(time);
        const amrex::Real ph   = wave.k * xv - wave.omega * time;
        const amrex::Real u_w  =  wave.a * r * wave.omega * std::sin(ph);
        const amrex::Real w_w  = -wave.a * r * wave.omega * std::cos(ph);
        const amrex::Real dhdy = amrex::Real(0.0);

        const amrex::Real vloc = vely_arr(i,j,k);
        const amrex::Real uloc = amrex::Real(0.25) * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                                     + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        const amrex::Real wloc = amrex::Real(0.25) * ( velz_arr(i,j  ,k) + velz_arr(i,j  ,k+1)
                                                     + velz_arr(i,j-1,k) + velz_arr(i,j-1,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;
        const amrex::Real w_r = wloc - w_w;

        const amrex::Real v_par = (v_r + w_r * dhdy) / std::sqrt(dhdy*dhdy + amrex::Real(1.0));
        const amrex::Real V0    = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);
        const amrex::Real rho   = amrex::Real(0.5) * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        if (stab) {
            const amrex::Real theta1 = amrex::Real(0.5) * ( cons_arr(i,j-1,k,RhoTheta_comp)/cons_arr(i,j-1,k,Rho_comp)
                                                          + cons_arr(i,j  ,k,RhoTheta_comp)/cons_arr(i,j  ,k,Rho_comp) );
            Cd = stab_cd(V0, theta1, sp);
        }
        return -rho * Cd * V0 * v_par;
    }

    // Surface heat flux rho<theta'w'> = -rho u_* theta_* at the cell center (i,j,k) at x = xc;
    // with a prescribed flux it is rho Q
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    amrex::Real
    compute_wave_t_flux (const int& i, const int& j, const int& k,
                         const amrex::Array4<const amrex::Real>& cons_arr,
                         const amrex::Array4<const amrex::Real>& velx_arr,
                         const amrex::Array4<const amrex::Real>& vely_arr,
                         const amrex::Array4<const amrex::Real>& velz_arr,
                         const amrex::Real& xc, const amrex::Real& time,
                         const MovingWave& wave, const StabParams& sp)
    {
        const amrex::Real r    = wave.ramp(time);
        const amrex::Real ph   = wave.k * xc - wave.omega * time;
        const amrex::Real u_w  =  wave.a * r * wave.omega * std::sin(ph);
        const amrex::Real w_w  = -wave.a * r * wave.omega * std::cos(ph);

        const amrex::Real uloc = amrex::Real(0.5) * ( velx_arr(i,j,k) + velx_arr(i+1,j,k) );
        const amrex::Real vloc = amrex::Real(0.5) * ( vely_arr(i,j,k) + vely_arr(i,j+1,k) );
        const amrex::Real wloc = amrex::Real(0.5) * ( velz_arr(i,j,k) + velz_arr(i,j,k+1) );

        const amrex::Real u_r = uloc - u_w;
        const amrex::Real v_r = vloc;
        const amrex::Real w_r = wloc - w_w;
        const amrex::Real V0  = std::sqrt(u_r*u_r + v_r*v_r + w_r*w_r) + amrex::Real(1.0e-15);

        const amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        const amrex::Real theta1 = cons_arr(i,j,k,RhoTheta_comp) / rho;

        amrex::Real ustar, tstar, den_m;
        iterate_ustar_tstar(V0, theta1, sp, ustar, tstar, den_m);

        return -rho * ustar * tstar;
    }

} // namespace WaveDrag
#endif

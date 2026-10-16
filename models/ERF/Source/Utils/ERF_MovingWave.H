#ifndef ERF_MOVING_WAVE_H_
#define ERF_MOVING_WAVE_H_

#include <cmath>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include "ERF_Constants.H"

/**
 * Monochromatic wave of the moving-terrain surfaces MovingSineWave, MovingSineWaveRest and
 * MovingSineWaveRamp, read from prob.Ampl, prob.wavelength and prob.wave_c.
 * The angular frequency is omega = c k for a prescribed phase speed c = prob.wave_c > 0,
 * else the deep-water value sqrt(g k).
 */
struct MovingWave
{
    amrex::Real a     = amrex::Real(0.0); // amplitude
    amrex::Real k     = amrex::Real(0.0); // wavenumber
    amrex::Real omega = amrex::Real(0.0); // angular frequency

    static MovingWave from_inputs ()
    {
        MovingWave wave;
        amrex::Real wavelength = amrex::Real(100.0);
        amrex::Real c          = amrex::Real(0.0);
        amrex::ParmParse pp("prob");
        pp.query("Ampl", wave.a);
        pp.query("wavelength", wavelength);
        pp.query("wave_c", c);
        wave.k     = amrex::Real(2.0) * PI / wavelength;
        wave.omega = (c > amrex::Real(0.0)) ? c * wave.k : std::sqrt(CONST_GRAV * wave.k);
        return wave;
    }

    // Ramp R(t) = (1 - cos(omega t))/2 for omega t < pi, 1 after, which starts the motion from rest
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real ramp (const amrex::Real t) const
    {
        return (omega * t < amrex::Real(PI)) ?
               amrex::Real(0.5) * (amrex::Real(1.0) - std::cos(omega * t)) : amrex::Real(1.0);
    }
};

#endif

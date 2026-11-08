#pragma once

#include <AMReX_REAL.H>

#include <algorithm>
#include <cmath>

#include "CRESCENT_Constants.H"

namespace crescent::init::stokes {

struct Wavevector
{
    int mode_x = 0;
    int mode_y = 0;
    amrex::Real kx = 0.0;
    amrex::Real ky = 0.0;
    amrex::Real magnitude = 0.0;
    amrex::Real actual_angle = 0.0;
};

inline void populatePhysicalComponents(
    Wavevector& wavevector,
    amrex::Real dk_x,
    amrex::Real dk_y)
{
    wavevector.kx = wavevector.mode_x * dk_x;
    wavevector.ky = wavevector.mode_y * dk_y;
    wavevector.magnitude = std::hypot(wavevector.kx, wavevector.ky);
    wavevector.actual_angle = std::atan2(wavevector.ky, wavevector.kx);
}

// Wavevector of the Stokes carrier on the periodic grid: |mode_x| = wave_mode_x
// away from the y axis, signs taken from the requested angle in every quadrant.
inline Wavevector computeWavevector(
    amrex::Real dk_x,
    amrex::Real dk_y,
    int wave_mode_x,
    amrex::Real wave_angle)
{
    constexpr amrex::Real axis_tolerance = 1.0e-10;
    const amrex::Real normalized_angle =
        std::remainder(wave_angle, crescent::constants::PI2);
    const amrex::Real cosine = std::cos(normalized_angle);
    const amrex::Real sine = std::sin(normalized_angle);

    Wavevector wavevector;
    if (std::abs(sine) < axis_tolerance) {
        wavevector.mode_x = cosine >= 0.0 ? wave_mode_x : -wave_mode_x;
    } else if (std::abs(cosine) < axis_tolerance) {
        wavevector.mode_y = sine >= 0.0 ? wave_mode_x : -wave_mode_x;
    } else {
        wavevector.mode_x = cosine > 0.0 ? wave_mode_x : -wave_mode_x;
        const amrex::Real target_ratio =
            cosine / sine * (dk_y / dk_x);
        wavevector.mode_y = static_cast<int>(
            std::round(wavevector.mode_x / target_ratio));
        if (wavevector.mode_y == 0) {
            wavevector.mode_y = sine > 0.0 ? 1 : -1;
        }
    }
    populatePhysicalComponents(wavevector, dk_x, dk_y);
    return wavevector;
}

} // namespace crescent::init::stokes

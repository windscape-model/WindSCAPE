#include <AMReX.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <iostream>

#include "CRESCENT_Constants.H"
#include "CRESCENT_StokesWavevector.H"
#include "CRESCENT_WaveState.H"

namespace {

bool hasModes(
    const crescent::init::stokes::Wavevector& wavevector,
    int mode_x,
    int mode_y)
{
    return wavevector.mode_x == mode_x && wavevector.mode_y == mode_y;
}

int runTest()
{
    using amrex::Real;
    using crescent::constants::PI;
    namespace stokes = crescent::init::stokes;

    constexpr Real dk = 0.1;
    constexpr int mode = 16;
    bool ok = true;

    // The x-mode magnitude is kept away from the y axis and the signs follow
    // the requested angle in every quadrant.
    const auto forward = stokes::computeWavevector(dk, dk, mode, PI / Real(6.0));
    const auto negative = stokes::computeWavevector(dk, dk, mode, -PI / Real(6.0));
    const auto backward =
        stokes::computeWavevector(dk, dk, mode, Real(5.0) * PI / Real(6.0));
    const auto minus_y = stokes::computeWavevector(dk, dk, mode, -PI / Real(2.0));
    const auto minus_x = stokes::computeWavevector(dk, dk, mode, PI);
    if (!hasModes(forward, 16, 9) ||
        !hasModes(negative, 16, -9) ||
        !hasModes(backward, -16, 9) ||
        !hasModes(minus_y, 0, -16) ||
        !hasModes(minus_x, -16, 0)) {
        std::cerr << "Stokes wavevector mapping does not cover all quadrants\n";
        ok = false;
    }

    // A positive crescent.peak_wavenumber_override replaces the initializer's
    // filter peak; the default -1 keeps it.
    crescent::shared::WaveRuntimeState wave;
    wave.spectral.peak_wavenumber_physical = 1.5;
    if (wave.spectralFilterSettings().peak_wavenumber_physical != 1.5) {
        std::cerr << "default filter peak is not the initializer's peak\n";
        ok = false;
    }
    wave.spectral.peak_wavenumber_override = 0.25;
    if (wave.spectralFilterSettings().peak_wavenumber_physical != 0.25) {
        std::cerr << "explicit peak-wavenumber override is inert\n";
        ok = false;
    }

    return ok ? 0 : 1;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    const int result = runTest();
    amrex::Finalize();
    return result;
}

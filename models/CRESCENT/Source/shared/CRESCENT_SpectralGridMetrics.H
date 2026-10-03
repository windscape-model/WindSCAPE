#ifndef CRESCENT_SPECTRAL_GRID_METRICS_H
#define CRESCENT_SPECTRAL_GRID_METRICS_H

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_Constants.H"
#include "CRESCENT_WaveState.H"

namespace crescent::shared {

struct SpectralGridMetrics
{
    int nx = 0;
    int ny = 0;
    int i0 = 0;
    int j0 = 0;
    amrex::Real dx = 0.0;
    amrex::Real dy = 0.0;
    amrex::Real Lx = 0.0;
    amrex::Real Ly = 0.0;
    amrex::Real dk_x = 0.0;
    amrex::Real dk_y = 0.0;
    amrex::Real k_nyquist_x = 0.0;
    amrex::Real k_nyquist_y = 0.0;

    [[nodiscard]] amrex::Real minNyquist() const
    {
        return (k_nyquist_x < k_nyquist_y) ? k_nyquist_x : k_nyquist_y;
    }

    [[nodiscard]] WavenumberSpacing wavenumberSpacing() const
    {
        return {dk_x, dk_y};
    }
};

inline SpectralGridMetrics makeSpectralGridMetrics(
    const amrex::Box& grid_box,
    amrex::Real dx,
    amrex::Real dy)
{
    SpectralGridMetrics metrics;
    metrics.nx = grid_box.length(0);
    metrics.ny = grid_box.length(1);
    metrics.i0 = grid_box.smallEnd(0);
    metrics.j0 = grid_box.smallEnd(1);
    metrics.dx = dx;
    metrics.dy = dy;
    metrics.Lx = metrics.nx * dx;
    metrics.Ly = metrics.ny * dy;
    metrics.dk_x = crescent::constants::PI2 / metrics.Lx;
    metrics.dk_y = crescent::constants::PI2 / metrics.Ly;
    metrics.k_nyquist_x = 0.5 * metrics.nx * metrics.dk_x;
    metrics.k_nyquist_y = 0.5 * metrics.ny * metrics.dk_y;
    return metrics;
}

inline SpectralGridMetrics makeSpectralGridMetrics(
    const amrex::Box& grid_box,
    const amrex::Geometry& geom)
{
    return makeSpectralGridMetrics(
        grid_box,
        geom.CellSize(0),
        geom.CellSize(1));
}

inline SpectralGridMetrics makeSpectralGridMetrics(
    const amrex::Box& grid_box,
    const amrex::GeometryData& geomdata)
{
    return makeSpectralGridMetrics(
        grid_box,
        geomdata.CellSize(0),
        geomdata.CellSize(1));
}

inline SpectralGridMetrics makeSpectralGridMetrics(
    const amrex::MultiFab& mf,
    const amrex::Geometry& geom)
{
    return makeSpectralGridMetrics(mf.boxArray().minimalBox(), geom);
}

inline SpectralGridMetrics makeSpectralGridMetrics(
    const amrex::MultiFab& mf,
    const amrex::GeometryData& geomdata)
{
    return makeSpectralGridMetrics(mf.boxArray().minimalBox(), geomdata);
}

} // namespace crescent::shared

#endif // CRESCENT_SPECTRAL_GRID_METRICS_H

#ifndef CRESCENT_SPECTRAL_WORKSPACE_H
#define CRESCENT_SPECTRAL_WORKSPACE_H

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "CRESCENT_SpectralWeightTable.H"

namespace crescent::shared {

struct SpectralWorkspace
{
    amrex::Vector<std::unique_ptr<amrex::FFT::R2C<amrex::Real>>> fft;

    amrex::Vector<amrex::cMultiFab> phi_expansion_spectrum;
    amrex::Vector<amrex::MultiFab> compute_rhs_metric;
    amrex::Vector<amrex::MultiFab> compute_rhs_phi_z_sq;
    amrex::Vector<amrex::MultiFab> expansion_accum;
    amrex::Vector<amrex::MultiFab> spectral_work;
    amrex::Vector<amrex::MultiFab> expansion_deriv;
    amrex::Vector<amrex::cMultiFab> deriv_spectral;

    // Pre-allocated scratch for directional derivative FFT pairs. Reused every
    // substage by executeDirectionalDerivativePair[PreservingGrowCells]; holds
    // the input field spectrum and the x/y derivative spectra.
    amrex::Vector<amrex::cMultiFab> directional_input_spectral;
    amrex::Vector<amrex::cMultiFab> directional_x_spectral;
    amrex::Vector<amrex::cMultiFab> directional_y_spectral;
    // Valid-region real-space scratch used by the growth-cell-preserving
    // directional derivative path so that FFT backward does not overwrite the
    // caller's ghost cells.
    amrex::Vector<amrex::MultiFab> directional_x_valid;
    amrex::Vector<amrex::MultiFab> directional_y_valid;
    // Pre-allocated scratch for crescent::spectral::dealiasing. Sized per level
    // to match the largest nGrow of dealiasing targets at that level
    // (scalar_nghost = 1 on L0, 2 on L1+). Zeroed before each use so that
    // ghost cells match the fresh-allocation semantics of the original
    // implementation.
    amrex::Vector<amrex::MultiFab> dealiasing_scratch;
    // |k|^n (n = 1..M) of the vertical-velocity series, and the level's
    // filter weights (parent low-pass, child combined filter).
    amrex::Vector<SpectralWeightTable> vertical_velocity_kpow;
    amrex::Vector<SpectralWeightTable> filter_weights;

    void resize(int nlevs_max)
    {
        fft.resize(nlevs_max);
        phi_expansion_spectrum.resize(nlevs_max);
        compute_rhs_metric.resize(nlevs_max);
        compute_rhs_phi_z_sq.resize(nlevs_max);
        expansion_accum.resize(nlevs_max);
        spectral_work.resize(nlevs_max);
        expansion_deriv.resize(nlevs_max);
        deriv_spectral.resize(nlevs_max);
        directional_input_spectral.resize(nlevs_max);
        directional_x_spectral.resize(nlevs_max);
        directional_y_spectral.resize(nlevs_max);
        directional_x_valid.resize(nlevs_max);
        directional_y_valid.resize(nlevs_max);
        dealiasing_scratch.resize(nlevs_max);
        vertical_velocity_kpow.resize(nlevs_max);
        filter_weights.resize(nlevs_max);
    }
};

} // namespace crescent::shared

#endif // CRESCENT_SPECTRAL_WORKSPACE_H

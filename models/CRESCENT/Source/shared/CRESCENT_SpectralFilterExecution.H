#pragma once

#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace crescent::spectral {

template <typename SpectralFilterOperation>
void applyInPlaceSpectralFilter(
    amrex::MultiFab& field,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    SpectralFilterOperation spectral_filter)
{
    fft->forwardThenBackward(field, field, spectral_filter);
}

template <typename SpectralFilterOperation>
void applySurfacePairSpectralFilter(
    amrex::MultiFab& eta,
    amrex::MultiFab& phi,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    SpectralFilterOperation spectral_filter)
{
    applyInPlaceSpectralFilter(eta, fft, spectral_filter);
    applyInPlaceSpectralFilter(phi, fft, spectral_filter);
}

} // namespace crescent::spectral

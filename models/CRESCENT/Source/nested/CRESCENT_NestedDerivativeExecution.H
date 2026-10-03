#ifndef CRESCENT_NESTED_DERIVATIVE_EXECUTION_H
#define CRESCENT_NESTED_DERIVATIVE_EXECUTION_H

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_FreeSurfaceViews.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_NestedSurfaceDerivatives.H"

namespace crescent::nested::derivative_execution {

struct PreparedPeriodicDerivativeContext
{
    shared::LevelWaveOperatorView level_wave_operator;
    const amrex::Geometry& level_geometry;
    shared::ConstFreeSurfaceStateView linear_derivative_state;
    shared::ConstFreeSurfaceStateView dno_surface_state;
    const amrex::Box& spectral_extract_box;
    int linear_derivative_spectral_filter_mode;
    int derivative_verbosity;
    int dno_verbosity;
    shared::DerivativeBundleView derivative_bundle;
    const amrex::MultiFab& wavenumber_pow_source;
    amrex::MultiFab& eta_pow_storage;
    amrex::MultiFab& phi_expansion_storage;
};

inline void computePreparedPeriodicDerivativeBundle(
    const PreparedPeriodicDerivativeContext& context)
{
    const int lev = context.level_wave_operator.lev;
    const auto& spectral_settings =
        context.level_wave_operator.spectral_settings;
    auto& spectral_workspace =
        context.level_wave_operator.spectral_workspace;

    crescent::nested::derivatives::computeSurfaceDerivatives(
        context.derivative_bundle.eta_x,
        context.derivative_bundle.eta_y,
        context.derivative_bundle.phi_x,
        context.derivative_bundle.phi_y,
        context.linear_derivative_state.eta,
        context.linear_derivative_state.phi,
        context.spectral_extract_box,
        context.level_geometry,
        context.level_wave_operator.fft,
        context.linear_derivative_spectral_filter_mode,
        spectral_settings.filter.peak_wavenumber_physical,
        spectral_settings.filter.beta1_filter,
        spectral_settings.filter.beta2_filter,
        spectral_workspace.directional_input_spectral[lev],
        spectral_workspace.directional_x_spectral[lev],
        spectral_workspace.directional_y_spectral[lev],
        spectral_workspace.directional_x_valid[lev],
        spectral_workspace.directional_y_valid[lev],
        context.derivative_verbosity);

    crescent::dno::computeNestedL1DNO(
        crescent::dno::makeSurfaceStateDNOContext(
            context.level_wave_operator,
            context.dno_verbosity,
            context.dno_surface_state,
            context.wavenumber_pow_source,
            context.eta_pow_storage,
            context.phi_expansion_storage,
            context.derivative_bundle.phi_z));
}

} // namespace crescent::nested::derivative_execution

#endif // CRESCENT_NESTED_DERIVATIVE_EXECUTION_H

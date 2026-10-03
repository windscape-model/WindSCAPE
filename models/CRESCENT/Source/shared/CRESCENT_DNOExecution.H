#ifndef CRESCENT_DNO_EXECUTION_H
#define CRESCENT_DNO_EXECUTION_H

#include <memory>

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_FreeSurfaceViews.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_SolverAccess.H"

namespace crescent::dno {

struct SurfaceStateDNOContext
{
    crescent::shared::LevelWaveOperatorView level_wave_operator;
    int dno_verbosity;
    crescent::shared::ConstFreeSurfaceStateView source_state;
    const amrex::MultiFab& wavenumber_pow_source;
    amrex::MultiFab& eta_pow_storage;
    amrex::MultiFab& phi_expansion_storage;
    amrex::MultiFab& phi_z_out;
};

SurfaceStateDNOContext makeSurfaceStateDNOContext(
    crescent::shared::LevelWaveOperatorView level_wave_operator,
    int dno_verbosity,
    crescent::shared::ConstFreeSurfaceStateView source_state,
    const amrex::MultiFab& wavenumber_pow_source,
    amrex::MultiFab& eta_pow_storage,
    amrex::MultiFab& phi_expansion_storage,
    amrex::MultiFab& phi_z_out);

void computeSolutionStateDNO(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const crescent::ownership::CoreSolverAccess& core,
    amrex::MultiFab& phi_z_out);

void computeSolutionStateDerivativeBundle(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const crescent::ownership::CoreSolverAccess& core,
    crescent::shared::DerivativeBundleView derivative_bundle);

void computeNestedL1DNO(const SurfaceStateDNOContext& context);

} // namespace crescent::dno

#endif // CRESCENT_DNO_EXECUTION_H

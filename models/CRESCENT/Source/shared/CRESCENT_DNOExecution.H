#ifndef CRESCENT_DNO_EXECUTION_H
#define CRESCENT_DNO_EXECUTION_H

#include <AMReX_MultiFab.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_LevelWaveOperatorView.H"

namespace crescent {
struct CoreSolverState;
}

namespace crescent::dno {

// Surface derivatives and the HOS vertical velocity of the level's solution
// state, on a periodic level whose transform box is the problem domain (the
// parent, or a single-level run).
void computeSolutionStateDerivativeBundle(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    crescent::CoreSolverState& core,
    crescent::shared::DerivativeBundleView derivative_bundle);

// HOS vertical velocity of a periodic child-level pair (eta, phi), whose
// transform box is the L1 box. The correction-field coupling projects the
// result itself (taper, project, taper), so the final projection is optional.
void computeNestedL1DNO(
    const crescent::shared::LevelWaveOperatorView& level_wave_operator,
    const amrex::MultiFab& eta,
    const amrex::MultiFab& phi,
    const amrex::MultiFab& wavenumber_pow,
    amrex::MultiFab& eta_pow,
    amrex::MultiFab& phi_expansion,
    amrex::MultiFab& phi_z,
    bool apply_final_projection);

} // namespace crescent::dno

#endif // CRESCENT_DNO_EXECUTION_H

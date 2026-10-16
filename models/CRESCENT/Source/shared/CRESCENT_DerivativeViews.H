#ifndef CRESCENT_DERIVATIVE_VIEWS_H
#define CRESCENT_DERIVATIVE_VIEWS_H

#include <AMReX_MultiFab.H>

namespace crescent::shared {

struct DerivativeBundleView {
    amrex::MultiFab& eta_x;
    amrex::MultiFab& eta_y;
    amrex::MultiFab& phi_x;
    amrex::MultiFab& phi_y;
    amrex::MultiFab& phi_z;
};

} // namespace crescent::shared

#endif // CRESCENT_DERIVATIVE_VIEWS_H

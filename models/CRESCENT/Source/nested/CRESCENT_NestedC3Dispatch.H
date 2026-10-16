#ifndef CRESCENT_NESTED_C3_DISPATCH_H
#define CRESCENT_NESTED_C3_DISPATCH_H

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "CRESCENT_NestedC3Periodicity.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"

namespace crescent::nested::c3 {

inline C3InterfaceRepairOptions makeInterfaceRepairOptions(
    int lev,
    const crescent::nested::AlgorithmConfig& nested_config,
    const crescent::nested::CouplingZoneState& coupling_zones)
{
    C3InterfaceRepairOptions repair_options;
    repair_options.buffer_width_x =
        coupling_zones.c3_buffer_width_x_computed[lev];
    repair_options.buffer_width_y =
        coupling_zones.c3_buffer_width_y_computed[lev];
    repair_options.boundary_penalty_weight = nested_config.c3_boundary_penalty;
    repair_options.coefficient_regularization_weight =
        nested_config.c3_coefficient_regularization;
    return repair_options;
}

inline void applyRepairPair(
    int lev,
    amrex::MultiFab& eta_mf,
    amrex::MultiFab& phi_mf,
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::BoxArray>& grids,
    const crescent::nested::AlgorithmConfig& nested_config,
    const crescent::nested::CouplingZoneState& coupling_zones,
    crescent::nested::C3RepairWorkspace& workspace)
{
    const auto repair_options =
        makeInterfaceRepairOptions(lev, nested_config, coupling_zones);

    crescent::nested::c3::applyC3InterfaceSeamRepair(
        eta_mf, geom[lev], grids[lev], repair_options, workspace);
    crescent::nested::c3::applyC3InterfaceSeamRepair(
        phi_mf, geom[lev], grids[lev], repair_options, workspace);
}

} // namespace crescent::nested::c3

#endif // CRESCENT_NESTED_C3_DISPATCH_H

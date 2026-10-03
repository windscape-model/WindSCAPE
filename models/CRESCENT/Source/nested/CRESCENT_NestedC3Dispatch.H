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
    const crescent::nested::CouplingZoneState& coupling_zones,
    int repair_verbosity)
{
    C3InterfaceRepairOptions repair_options;
    repair_options.buffer_width = coupling_zones.c3_buffer_width_computed[lev];
    repair_options.verbosity = repair_verbosity;
    repair_options.apply_x = (nested_config.c3_apply_x == 1);
    repair_options.apply_y = (nested_config.c3_apply_y == 1);
    repair_options.boundary_penalty_weight = nested_config.c3_boundary_penalty;
    repair_options.coefficient_regularization_weight =
        nested_config.c3_coefficient_regularization;
    repair_options.x_pass_y_smoothing_passes =
        nested_config.c3_x_pass_y_smoothing_passes;
    repair_options.x_pass_corner_taper = nested_config.c3_x_pass_corner_taper;
    repair_options.x_pass_alpha_y_smoothing_passes =
        nested_config.c3_x_alpha_y_smoothing_passes;
    repair_options.x_alpha_y_regularization_weight =
        nested_config.c3_x_alpha_y_regularization;
    repair_options.y_alpha_x_regularization_weight =
        nested_config.c3_y_alpha_x_regularization;
    repair_options.x_alpha_y_mean_regularization_weight =
        nested_config.c3_x_alpha_y_mean_regularization;
    repair_options.y_alpha_x_mean_regularization_weight =
        nested_config.c3_y_alpha_x_mean_regularization;
    repair_options.x_transverse_harmonics = nested_config.c3_x_transverse_harmonics;
    repair_options.y_transverse_harmonics = nested_config.c3_y_transverse_harmonics;
    // Preserve current production behavior by leaving projection harmonics at
    // their default values in this paired repair path.
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
    int repair_verbosity,
    crescent::nested::C3RepairWorkspace& workspace)
{
    const auto repair_options =
        makeInterfaceRepairOptions(lev, nested_config, coupling_zones, repair_verbosity);

    crescent::nested::c3::applyC3InterfaceSeamRepair(
        eta_mf, geom[lev], grids[lev], repair_options, workspace);
    crescent::nested::c3::applyC3InterfaceSeamRepair(
        phi_mf, geom[lev], grids[lev], repair_options, workspace);
}

} // namespace crescent::nested::c3

#endif // CRESCENT_NESTED_C3_DISPATCH_H

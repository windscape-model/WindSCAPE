#include "CRESCENT_NestedAlgorithmInputs.H"

#include "CRESCENT_RetiredOptions.H"
#include "CRESCENT_WaveState.H"

#include <cmath>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

namespace {

void readRefinementRegionParameters(
    crescent::nested::AlgorithmConfig& nested_config,
    crescent::nested::config::UserSpecifiedOptions& user_specified_options,
    amrex::ParmParse& pp)
{
    const bool has_l1_box_lo = pp.contains("L1_box_lo");
    const bool has_l1_box_hi = pp.contains("L1_box_hi");
    user_specified_options.l1_refinement_box = has_l1_box_lo && has_l1_box_hi;

    if (!user_specified_options.l1_refinement_box) {
        return;
    }

    amrex::Vector<amrex::Real> l1_box_lo(2);
    amrex::Vector<amrex::Real> l1_box_hi(2);
    pp.getarr("L1_box_lo", l1_box_lo, 0, 2);
    pp.getarr("L1_box_hi", l1_box_hi, 0, 2);

    nested_config.l1_box_x_lo = l1_box_lo[0];
    nested_config.l1_box_y_lo = l1_box_lo[1];
    nested_config.l1_box_x_hi = l1_box_hi[0];
    nested_config.l1_box_y_hi = l1_box_hi[1];
}

using crescent::shared::requireRetiredValue;

void readBlendingParameters(amrex::ParmParse& pp)
{
    // Only the septic tensor-product mask remains (former mode 3, alpha 1).
    requireRetiredValue(pp, "blending_function_mode", 3);
    requireRetiredValue(pp, "use_sequential_blending", 0);
    requireRetiredValue(pp, "correct_sequential_blending", 1);
    requireRetiredValue(pp, "blending_tensor_alpha", amrex::Real(1.0));
    requireRetiredValue(pp, "blending_superellipse_gamma", amrex::Real(1.0));
    requireRetiredValue(pp, "blending_superellipse_p", amrex::Real(4.0));
    requireRetiredValue(pp, "damping_steepness", amrex::Real(2.0));
    requireRetiredValue(pp, "blending_correction_y_smoothing_passes", 0);
}

void readC3Parameters(
    crescent::nested::AlgorithmConfig& nested_config,
    amrex::ParmParse& pp)
{
    for (const char* name : {"use_c3_repair", "c3_use_interface", "c3_apply_x",
                             "c3_apply_y", "c3_directional_buffer_widths"}) {
        requireRetiredValue(pp, name, 1);
    }
    for (const char* name : {"c3_x_pass_y_smoothing_passes", "c3_x_pass_corner_taper",
                             "c3_x_alpha_y_smoothing_passes", "c3_repair_after_filter",
                             "c3_x_transverse_harmonics", "c3_y_transverse_harmonics",
                             "c3_x_alpha_y_projection_harmonics",
                             "c3_y_alpha_x_projection_harmonics"}) {
        requireRetiredValue(pp, name, 0);
    }
    for (const char* name : {"c3_x_alpha_y_regularization", "c3_y_alpha_x_regularization",
                             "c3_x_alpha_y_mean_regularization",
                             "c3_y_alpha_x_mean_regularization"}) {
        requireRetiredValue(pp, name, amrex::Real(0.0));
    }
    // The buffer width is the distance from the child box edge to Omega_S.
    requireRetiredValue(pp, "c3_buffer_width", 24);
    pp.query("c3_boundary_penalty", nested_config.c3_boundary_penalty);
    pp.query("c3_coefficient_regularization", nested_config.c3_coefficient_regularization);
}

} // namespace

namespace crescent::nested::config {

ParsedAlgorithmInputs readAlgorithmInputs(amrex::ParmParse& pp)
{
    ParsedAlgorithmInputs inputs;
    readRefinementRegionParameters(
        inputs.config, inputs.user_specified_options, pp);

    pp.query("use_crescent_algorithm", inputs.config.use_crescent_algorithm);

    inputs.user_specified_options.blending_steepness =
        pp.contains("blending_steepness");

    pp.query("roi_accurate_fraction", inputs.config.roi_accurate_fraction);
    pp.query("roi_spectral_fraction", inputs.config.roi_spectral_fraction);
    pp.query("state_interpolation_order", inputs.config.state_interpolation_order);
    pp.query("derivative_interpolation_order", inputs.config.derivative_interpolation_order);
    pp.query("l1_coupling_mode", inputs.config.l1_coupling_mode);
    requireRetiredValue(pp, "l1_rhs_eta_source", 0);
    requireRetiredValue(pp, "ccf_band_restrict", 1);
    requireRetiredValue(pp, "enable_final_output_blended_derivative_filter", 0);

    for (const char* name : {"ccf_background_tendency", "ccf_linear_consistency",
                             "ccf_adjustment_consistency", "ccf_transform_admissibility",
                             "ccf_w1_projection_consistency", "l1_final_dno_projection",
                             "l1_correct_bilinear_interpolation",
                             "filter_periodic_copies_before_linear_derivatives",
                             "l1_filter_preserve_mean", "consistent_applied_pressure",
                             "enable_post_rk4_state_filter_and_blend",
                             "consistent_output_diagnostics"}) {
        requireRetiredValue(pp, name, 1);
    }
    for (const char* name : {"linear_derivative_spectral_filter_mode",
                             "l1_stage_phi_z_cleanup",
                             "enable_post_timestep_low_pass_filter"}) {
        requireRetiredValue(pp, name, 0);
    }
    readBlendingParameters(pp);
    readC3Parameters(inputs.config, pp);

    return inputs;
}

void validateAlgorithmSelection(
    const AlgorithmConfig& nested_config,
    int max_level)
{
    if (nested_config.use_crescent_algorithm != 0 &&
        nested_config.use_crescent_algorithm != 1) {
        amrex::Abort("crescent.use_crescent_algorithm must be 0 or 1");
    }
    if (max_level > ActiveNestedL1Level) {
        amrex::Abort(
            "CRESCENT supports at most one nested level; amr.max_level must be 0 or 1");
    }
    if (max_level >= 1 && nested_config.use_crescent_algorithm != 1) {
        amrex::Abort(
            "Nested L1 runs require crescent.use_crescent_algorithm = 1");
    }
    if (max_level == 0 && nested_config.use_crescent_algorithm != 0) {
        amrex::Abort(
            "crescent.use_crescent_algorithm=1 requires amr.max_level=1");
    }
    if (max_level >= 1 && amrex::ParallelDescriptor::NProcs() != 1) {
        amrex::Abort("Nested CRESCENT runs currently require exactly one MPI rank");
    }
}

void validateRefinementRegionInputs(
    const AlgorithmConfig& nested_config,
    const UserSpecifiedOptions& user_specified_options,
    int max_level)
{
    if (max_level < 1) {
        return;
    }

    if (!user_specified_options.l1_refinement_box) {
        amrex::Abort(
            "Nested L1 runs require both crescent.L1_box_lo and crescent.L1_box_hi");
    }

    if (!std::isfinite(nested_config.l1_box_x_lo) ||
        !std::isfinite(nested_config.l1_box_x_hi) ||
        !std::isfinite(nested_config.l1_box_y_lo) ||
        !std::isfinite(nested_config.l1_box_y_hi) ||
        !(nested_config.l1_box_x_lo < nested_config.l1_box_x_hi) ||
        !(nested_config.l1_box_y_lo < nested_config.l1_box_y_hi)) {
        amrex::Abort(
            "Nested L1 refinement box must satisfy "
            "L1_box_lo[x,y] < L1_box_hi[x,y]");
    }
}

void validateInterpolationInputs(const AlgorithmConfig& nested_config)
{
    if (nested_config.state_interpolation_order != 1 &&
        nested_config.state_interpolation_order != 2) {
        amrex::Abort(
            "crescent.state_interpolation_order must be 1 (bilinear) or 2 (quadratic)");
    }
    if (nested_config.derivative_interpolation_order != 1 &&
        nested_config.derivative_interpolation_order != 2) {
        amrex::Abort(
            "crescent.derivative_interpolation_order must be 1 (bilinear) or 2 (quadratic)");
    }
}

void validateBlendingInputs(const AlgorithmConfig& nested_config)
{
    if (!std::isfinite(nested_config.roi_spectral_fraction) ||
        !std::isfinite(nested_config.roi_accurate_fraction) ||
        !(nested_config.roi_spectral_fraction > 0.0) ||
        !(nested_config.roi_spectral_fraction <
          nested_config.roi_accurate_fraction) ||
        !(nested_config.roi_accurate_fraction < 0.5)) {
        amrex::Abort(
            "CRESCENT ROI fractions must satisfy 0 < roi_spectral_fraction "
            "< roi_accurate_fraction < 0.5");
    }
}

void validateC3Inputs(const AlgorithmConfig& nested_config)
{
    if (nested_config.c3_boundary_penalty < 0.0 ||
        !std::isfinite(nested_config.c3_boundary_penalty)) {
        amrex::Abort(
            "crescent.c3_boundary_penalty must be non-negative and finite");
    }
    if (nested_config.c3_coefficient_regularization < 0.0 ||
        !std::isfinite(nested_config.c3_coefficient_regularization)) {
        amrex::Abort(
            "crescent.c3_coefficient_regularization must be non-negative and finite");
    }
    if (nested_config.l1_coupling_mode != 1 &&
        nested_config.l1_coupling_mode != 2) {
        amrex::Abort(
            "crescent.l1_coupling_mode must be 1 (total field) or 2 "
            "(correction field); the derivative blend (0) was retired and is "
            "available at commit f39cacf and earlier");
    }
    if (nested_config.l1_coupling_mode == 2) {
        // Correction-field requirements (docs/R2_CORRECTION_FIELD_SPEC.md §6.5).
        if (nested_config.state_interpolation_order !=
            nested_config.derivative_interpolation_order) {
            amrex::Abort(
                "l1_coupling_mode=2 requires state_interpolation_order == "
                "derivative_interpolation_order (single interpolation operator)");
        }
    }
}

void emitInputWarnings(
    const AlgorithmConfig& nested_config,
    const UserSpecifiedOptions& user_specified_options)
{
    if (nested_config.use_crescent_algorithm == 1 &&
        user_specified_options.blending_steepness) {
        amrex::Print()
            << "WARNING: crescent.blending_steepness has no effect: the "
               "blending mask uses the fixed septic smoothstep."
            << std::endl;
    }
}

} // namespace crescent::nested::config

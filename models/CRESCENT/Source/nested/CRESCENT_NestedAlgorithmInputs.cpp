#include "CRESCENT_NestedAlgorithmInputs.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

namespace {

void emitBlendingSteepnessWarning(
    const crescent::nested::AlgorithmConfig& nested_config)
{
    const int active_blending_mode =
        crescent::nested::resolveActiveBlendingMode(nested_config);

    amrex::Print()
        << "WARNING: crescent.blending_steepness is accepted but ignored by the active CRESCENT sigma builder "
        << "(mode " << active_blending_mode << ": "
        << crescent::nested::blendingModeName(active_blending_mode) << "); ";

    switch (active_blending_mode) {
    case 0:
    case 1:
        amrex::Print()
            << "this mode uses a fixed septic smoothstep transition between Omega_A and Omega_S.";
        break;
    case 2:
        amrex::Print()
            << "this mode uses a harmonic sigma solve with fixed Omega_A/Omega_S endpoint constraints.";
        break;
    case 3:
        amrex::Print()
            << "shape control comes from crescent.blending_tensor_alpha.";
        break;
    case 4:
        amrex::Print()
            << "this mode uses a harmonic sigma solve followed by a fixed septic smoothstep remap.";
        break;
    case 5:
        amrex::Print()
            << "shape control comes from crescent.blending_superellipse_gamma and crescent.blending_superellipse_p.";
        break;
    default:
        amrex::Print()
            << "no active sigma builder currently uses this knob.";
        break;
    }

    amrex::Print() << std::endl;
}

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

void readBlendingParameters(
    crescent::nested::AlgorithmConfig& nested_config,
    amrex::ParmParse& pp)
{
    pp.query("blending_tensor_alpha", nested_config.blending_tensor_alpha);
    pp.query("blending_superellipse_gamma", nested_config.blending_superellipse_gamma);
    pp.query("blending_superellipse_p", nested_config.blending_superellipse_p);
    pp.query(
        "blending_correction_y_smoothing_passes",
        nested_config.blending_correction_y_smoothing_passes);
}

void readC3Parameters(
    crescent::nested::AlgorithmConfig& nested_config,
    amrex::ParmParse& pp)
{
    pp.query("use_c3_repair", nested_config.use_c3_repair);
    pp.query("c3_apply_x", nested_config.c3_apply_x);
    pp.query("c3_apply_y", nested_config.c3_apply_y);
    pp.query("c3_x_pass_y_smoothing_passes", nested_config.c3_x_pass_y_smoothing_passes);
    pp.query("c3_x_pass_corner_taper", nested_config.c3_x_pass_corner_taper);
    pp.query("c3_x_alpha_y_smoothing_passes", nested_config.c3_x_alpha_y_smoothing_passes);
    pp.query("c3_x_alpha_y_regularization", nested_config.c3_x_alpha_y_regularization);
    pp.query("c3_y_alpha_x_regularization", nested_config.c3_y_alpha_x_regularization);
    pp.query(
        "c3_x_alpha_y_mean_regularization",
        nested_config.c3_x_alpha_y_mean_regularization);
    pp.query(
        "c3_y_alpha_x_mean_regularization",
        nested_config.c3_y_alpha_x_mean_regularization);
    pp.query("c3_repair_after_filter", nested_config.c3_repair_after_filter);
    pp.query("c3_x_transverse_harmonics", nested_config.c3_x_transverse_harmonics);
    pp.query("c3_y_transverse_harmonics", nested_config.c3_y_transverse_harmonics);
    pp.query(
        "c3_x_alpha_y_projection_harmonics",
        nested_config.c3_x_alpha_y_projection_harmonics);
    pp.query(
        "c3_y_alpha_x_projection_harmonics",
        nested_config.c3_y_alpha_x_projection_harmonics);
    pp.query("c3_use_interface", nested_config.c3_use_interface);
    pp.query("c3_buffer_width", nested_config.c3_buffer_width);
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

    inputs.user_specified_options.damping_steepness =
        pp.contains("damping_steepness");
    inputs.user_specified_options.blending_steepness =
        pp.contains("blending_steepness");

    pp.query("roi_accurate_fraction", inputs.config.roi_accurate_fraction);
    pp.query("roi_spectral_fraction", inputs.config.roi_spectral_fraction);
    pp.query("damping_steepness", inputs.config.damping_steepness);
    pp.query("blending_steepness", inputs.config.blending_steepness);
    pp.query("state_interpolation_order", inputs.config.state_interpolation_order);
    pp.query("derivative_interpolation_order", inputs.config.derivative_interpolation_order);
    pp.query("blending_function_mode", inputs.config.blending_function_mode);
    pp.query("use_sequential_blending", inputs.config.use_sequential_blending);
    pp.query(
        "filter_periodic_copies_before_linear_derivatives",
        inputs.config.filter_periodic_copies_before_linear_derivatives);
    pp.query(
        "linear_derivative_spectral_filter_mode",
        inputs.config.linear_derivative_spectral_filter_mode);
    pp.query(
        "enable_post_timestep_low_pass_filter",
        inputs.config.enable_post_timestep_low_pass_filter);
    pp.query(
        "enable_post_rk4_state_filter_and_blend",
        inputs.config.enable_post_rk4_state_filter_and_blend);
    pp.query(
        "enable_final_output_blended_derivative_filter",
        inputs.config.enable_final_output_blended_derivative_filter);
    pp.query(
        "enable_initial_c3_effect_diagnostic",
        inputs.config.enable_initial_c3_effect_diagnostic);
    pp.query(
        "enable_nested_initialization_validation",
        inputs.config.enable_nested_initialization_validation);

    readBlendingParameters(inputs.config, pp);
    readC3Parameters(inputs.config, pp);

    return inputs;
}

void validateAlgorithmSelection(
    const AlgorithmConfig& nested_config,
    int max_level)
{
    if (max_level >= 1 && nested_config.use_crescent_algorithm != 1) {
        amrex::Abort(
            "Nested L1 runs require crescent.use_crescent_algorithm = 1");
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

    if (!(nested_config.l1_box_x_lo < nested_config.l1_box_x_hi) ||
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
    if (nested_config.blending_function_mode < 0 ||
        nested_config.blending_function_mode > 5) {
        amrex::Abort(
            "crescent.blending_function_mode must be 0 (Euclidean), 1 (sequential), "
            "2 (harmonic), 3 (tensor-product), 4 (harmonic+septic), or 5 (superellipse)");
    }
    if (nested_config.use_sequential_blending != 0 &&
        nested_config.use_sequential_blending != 1) {
        amrex::Abort("crescent.use_sequential_blending must be 0 or 1");
    }
    if (nested_config.blending_tensor_alpha <= 0.0) {
        amrex::Abort("crescent.blending_tensor_alpha must be > 0");
    }
    if (nested_config.blending_superellipse_gamma <= 0.0) {
        amrex::Abort("crescent.blending_superellipse_gamma must be > 0");
    }
    if (nested_config.blending_superellipse_p <= 0.0) {
        amrex::Abort("crescent.blending_superellipse_p must be > 0");
    }
    if (nested_config.blending_correction_y_smoothing_passes < 0) {
        amrex::Abort("crescent.blending_correction_y_smoothing_passes must be >= 0");
    }
}

void validateC3Inputs(const AlgorithmConfig& nested_config)
{
    if (nested_config.c3_apply_x != 0 && nested_config.c3_apply_x != 1) {
        amrex::Abort("crescent.c3_apply_x must be 0 or 1");
    }
    if (nested_config.c3_apply_y != 0 && nested_config.c3_apply_y != 1) {
        amrex::Abort("crescent.c3_apply_y must be 0 or 1");
    }
    if (nested_config.filter_periodic_copies_before_linear_derivatives != 0 &&
        nested_config.filter_periodic_copies_before_linear_derivatives != 1) {
        amrex::Abort(
            "crescent.filter_periodic_copies_before_linear_derivatives must be 0 or 1");
    }
    if (nested_config.linear_derivative_spectral_filter_mode < 0 ||
        nested_config.linear_derivative_spectral_filter_mode > 2) {
        amrex::Abort("crescent.linear_derivative_spectral_filter_mode must be 0, 1, or 2");
    }
    if (nested_config.enable_post_timestep_low_pass_filter != 0 &&
        nested_config.enable_post_timestep_low_pass_filter != 1) {
        amrex::Abort("crescent.enable_post_timestep_low_pass_filter must be 0 or 1");
    }
    if (nested_config.enable_post_rk4_state_filter_and_blend != 0 &&
        nested_config.enable_post_rk4_state_filter_and_blend != 1) {
        amrex::Abort("crescent.enable_post_rk4_state_filter_and_blend must be 0 or 1");
    }
    if (nested_config.enable_final_output_blended_derivative_filter != 0 &&
        nested_config.enable_final_output_blended_derivative_filter != 1) {
        amrex::Abort(
            "crescent.enable_final_output_blended_derivative_filter must be 0 or 1");
    }
    if (nested_config.use_crescent_algorithm == 1 &&
        nested_config.use_c3_repair == 1 &&
        nested_config.c3_apply_x == 0 &&
        nested_config.c3_apply_y == 0) {
        amrex::Abort(
            "crescent.c3_apply_x and crescent.c3_apply_y cannot both be 0 "
            "when C3 repair is enabled");
    }
    if (nested_config.c3_x_pass_y_smoothing_passes < 0) {
        amrex::Abort("crescent.c3_x_pass_y_smoothing_passes must be >= 0");
    }
    if (nested_config.c3_x_pass_corner_taper != 0 &&
        nested_config.c3_x_pass_corner_taper != 1) {
        amrex::Abort("crescent.c3_x_pass_corner_taper must be 0 or 1");
    }
    if (nested_config.c3_x_alpha_y_smoothing_passes < 0) {
        amrex::Abort("crescent.c3_x_alpha_y_smoothing_passes must be >= 0");
    }
    if (nested_config.c3_x_alpha_y_regularization < 0.0) {
        amrex::Abort("crescent.c3_x_alpha_y_regularization must be >= 0");
    }
    if (nested_config.c3_y_alpha_x_regularization < 0.0) {
        amrex::Abort("crescent.c3_y_alpha_x_regularization must be >= 0");
    }
    if (nested_config.c3_x_alpha_y_mean_regularization < 0.0) {
        amrex::Abort("crescent.c3_x_alpha_y_mean_regularization must be >= 0");
    }
    if (nested_config.c3_y_alpha_x_mean_regularization < 0.0) {
        amrex::Abort("crescent.c3_y_alpha_x_mean_regularization must be >= 0");
    }
    if (nested_config.c3_repair_after_filter != 0 &&
        nested_config.c3_repair_after_filter != 1) {
        amrex::Abort("crescent.c3_repair_after_filter must be 0 or 1");
    }
    if (nested_config.c3_x_transverse_harmonics < 0) {
        amrex::Abort("crescent.c3_x_transverse_harmonics must be >= 0");
    }
    if (nested_config.c3_y_transverse_harmonics < 0) {
        amrex::Abort("crescent.c3_y_transverse_harmonics must be >= 0");
    }
    if (nested_config.c3_x_alpha_y_projection_harmonics < 0) {
        amrex::Abort("crescent.c3_x_alpha_y_projection_harmonics must be >= 0");
    }
    if (nested_config.c3_y_alpha_x_projection_harmonics < 0) {
        amrex::Abort("crescent.c3_y_alpha_x_projection_harmonics must be >= 0");
    }
}

void emitInputWarnings(
    const AlgorithmConfig& nested_config,
    const UserSpecifiedOptions& user_specified_options)
{
    if (nested_config.use_crescent_algorithm == 1 &&
        user_specified_options.blending_steepness) {
        emitBlendingSteepnessWarning(nested_config);
    }
    if (nested_config.use_crescent_algorithm == 1 &&
        nested_config.use_c3_repair == 1 &&
        user_specified_options.damping_steepness) {
        amrex::Print()
            << "WARNING: crescent.damping_steepness is ignored when CRESCENT C3 repair is active." << std::endl;
    }
}

} // namespace crescent::nested::config

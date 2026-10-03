#ifndef CRESCENT_NESTED_DERIVATIVE_COUPLING_H
#define CRESCENT_NESTED_DERIVATIVE_COUPLING_H

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_NestedDerivativeBundles.H"
#include "CRESCENT_NestedFilter.H"
#include "CRESCENT_Spectral.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"
#include "shared/CRESCENT_LevelWaveOperatorView.H"

namespace crescent::nested::derivative_coupling {

using DerivativeBundleView = crescent::shared::DerivativeBundleView;
using LevelWaveSpectralSettings = crescent::shared::LevelWaveSpectralSettings;
using LevelEvolutionContext = crescent::shared::LevelEvolutionContext;
using amrex::Geometry;
using amrex::MultiFab;
using amrex::Real;

struct NestedFinalOutputBlendContext
{
    int derivative_interpolation_order;
    int blending_correction_y_smoothing_passes;
    int enable_final_output_blended_derivative_filter;
    LevelWaveSpectralSettings spectral_settings;

    const Geometry& core_geom;
    const Geometry& nested_geom;
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft;
    const MultiFab& sigma_blend;
    DerivativeBundleView native_l1_bundle;
    MultiFab& core_stage_derivative_scratch;
    MultiFab& core_interp_patch_scratch;
    MultiFab& correction_delta_raw;
    MultiFab& correction_delta_smooth;
    MultiFab& phi_z_filter_scratch;
    MultiFab& dealiasing_scratch;
};

inline NestedFinalOutputBlendContext makeNestedFinalOutputBlendContext(
    int derivative_interpolation_order,
    int blending_correction_y_smoothing_passes,
    int enable_final_output_blended_derivative_filter,
    LevelWaveSpectralSettings spectral_settings,
    const Geometry& core_geom,
    const Geometry& nested_geom,
    std::unique_ptr<amrex::FFT::R2C<amrex::Real>>& fft,
    const MultiFab& sigma_blend,
    DerivativeBundleView native_l1_bundle,
    MultiFab& core_stage_derivative_scratch,
    MultiFab& core_interp_patch_scratch,
    MultiFab& correction_delta_raw,
    MultiFab& correction_delta_smooth,
    MultiFab& phi_z_filter_scratch,
    MultiFab& dealiasing_scratch)
{
    return {
        derivative_interpolation_order,
        blending_correction_y_smoothing_passes,
        enable_final_output_blended_derivative_filter,
        spectral_settings,
        core_geom,
        nested_geom,
        fft,
        sigma_blend,
        native_l1_bundle,
        core_stage_derivative_scratch,
        core_interp_patch_scratch,
        correction_delta_raw,
        correction_delta_smooth,
        phi_z_filter_scratch,
        dealiasing_scratch
    };
}

inline NestedFinalOutputBlendContext makeNestedFinalOutputBlendContext(
    LevelEvolutionContext& level_context,
    DerivativeBundleView native_l1_bundle)
{
    const int lev = level_context.level();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            level_context.coreAccess(),
            level_context.levelMesh().domain());
    const auto& nested_config = level_context.nestedConfig();
    auto& coupling_derivatives = level_context.nestedRuntime().coupling_derivatives;
    auto& periodic_derivative_scratch =
        level_context.nestedRuntime().periodic_derivative_scratch[lev];

    return makeNestedFinalOutputBlendContext(
        nested_config.derivative_interpolation_order,
        nested_config.blending_correction_y_smoothing_passes,
        nested_config.enable_final_output_blended_derivative_filter,
        level_wave_operator.spectral_settings,
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_wave_operator.fft,
        level_context.nestedRuntime().coupling_zones.sigma_blend[lev],
        native_l1_bundle,
        coupling_derivatives.core_stage_derivative_scratch,
        coupling_derivatives.core_interp_patch_scratch,
        coupling_derivatives.blend_correction_delta_raw,
        coupling_derivatives.blend_correction_delta_smooth,
        periodic_derivative_scratch.phi_z_filter_scratch,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev]);
}

inline void interpolateAndBlendCoreDerivativesToL1(
    DerivativeBundleView core_n_bundle,
    DerivativeBundleView core_np1_bundle,
    Real temporal_weight,
    const Geometry& core_geom,
    const Geometry& nested_geom,
    int derivative_interpolation_order,
    DerivativeBundleView native_l1_bundle,
    const MultiFab& sigma_blend,
    int blending_correction_y_smoothing_passes,
    MultiFab& core_stage_derivative_scratch,
    MultiFab& core_interp_patch_scratch,
    MultiFab& correction_delta_raw,
    MultiFab& correction_delta_smooth,
    DerivativeBundleView core_derivatives_on_l1_bundle,
    DerivativeBundleView blended_derivative_bundle)
{
    crescent::nested::bundles::interpolateCoreDerivativesToL1(
        core_n_bundle,
        core_np1_bundle,
        temporal_weight,
        core_geom,
        nested_geom,
        derivative_interpolation_order,
        core_stage_derivative_scratch,
        core_interp_patch_scratch,
        core_derivatives_on_l1_bundle);

    crescent::nested::bundles::blendL1DerivativesWithCore(
        native_l1_bundle,
        core_derivatives_on_l1_bundle,
        sigma_blend,
        blending_correction_y_smoothing_passes,
        correction_delta_raw,
        correction_delta_smooth,
        blended_derivative_bundle);
}

inline void filterNestedL1FinalOutputDerivatives(
    const NestedFinalOutputBlendContext& blend_context,
    DerivativeBundleView final_blended_bundle)
{
    if (blend_context.enable_final_output_blended_derivative_filter != 1) {
        return;
    }

    crescent::nested::filter::applyHouLiFilter_L1(
        final_blended_bundle.eta_x, final_blended_bundle.eta_y,
        blend_context.fft, blend_context.nested_geom.data(),
        blend_context.spectral_settings,
        false
    );
    crescent::nested::filter::applyHouLiFilter_L1(
        final_blended_bundle.phi_x, final_blended_bundle.phi_y,
        blend_context.fft, blend_context.nested_geom.data(),
        blend_context.spectral_settings,
        false
    );

    crescent::nested::filter::applyHouLiFilter_L1(
        final_blended_bundle.phi_z,
        blend_context.phi_z_filter_scratch,
        blend_context.fft,
        blend_context.nested_geom.data(),
        blend_context.spectral_settings,
        false
    );

    crescent::spectral::dealiasing(
        final_blended_bundle.eta_x,
        blend_context.fft,
        blend_context.nested_geom.data(),
        blend_context.dealiasing_scratch);
    crescent::spectral::dealiasing(
        final_blended_bundle.eta_y,
        blend_context.fft,
        blend_context.nested_geom.data(),
        blend_context.dealiasing_scratch);
    crescent::spectral::dealiasing(
        final_blended_bundle.phi_x,
        blend_context.fft,
        blend_context.nested_geom.data(),
        blend_context.dealiasing_scratch);
    crescent::spectral::dealiasing(
        final_blended_bundle.phi_y,
        blend_context.fft,
        blend_context.nested_geom.data(),
        blend_context.dealiasing_scratch);
    crescent::spectral::dealiasing(
        final_blended_bundle.phi_z,
        blend_context.fft,
        blend_context.nested_geom.data(),
        blend_context.dealiasing_scratch);
}

inline void prepareNestedL1FinalOutputDerivatives(
    const NestedFinalOutputBlendContext& blend_context,
    DerivativeBundleView core_final_bundle,
    DerivativeBundleView core_final_on_l1_bundle,
    DerivativeBundleView final_blended_bundle)
{
    interpolateAndBlendCoreDerivativesToL1(
        core_final_bundle,
        core_final_bundle,
        0.0,
        blend_context.core_geom,
        blend_context.nested_geom,
        blend_context.derivative_interpolation_order,
        blend_context.native_l1_bundle,
        blend_context.sigma_blend,
        blend_context.blending_correction_y_smoothing_passes,
        blend_context.core_stage_derivative_scratch,
        blend_context.core_interp_patch_scratch,
        blend_context.correction_delta_raw,
        blend_context.correction_delta_smooth,
        core_final_on_l1_bundle,
        final_blended_bundle);

    filterNestedL1FinalOutputDerivatives(
        blend_context,
        final_blended_bundle);
}

} // namespace crescent::nested::derivative_coupling

#endif // CRESCENT_NESTED_DERIVATIVE_COUPLING_H

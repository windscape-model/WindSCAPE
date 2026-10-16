#ifndef CRESCENT_NESTED_DERIVATIVE_COUPLING_H
#define CRESCENT_NESTED_DERIVATIVE_COUPLING_H

#include "CRESCENT_DerivativeViews.H"
#include "CRESCENT_NestedDerivativeBundles.H"
#include "CRESCENT_SolverContext.H"

namespace crescent::nested::derivative_coupling {

// Output bundle of the total-field child: the native periodic derivative
// bundle blended with the parent step-end bundle prolonged to the child.
// Every operation here is pointwise; no non-periodic field is transformed.
inline void prepareNestedL1FinalOutputDerivatives(
    crescent::shared::LevelContext& level_context,
    crescent::shared::DerivativeBundleView native_l1_bundle)
{
    const int lev = level_context.level();
    const auto& nested_config = level_context.nestedConfig();
    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;
    auto parent_on_l1 = coupling_derivatives.core_derivatives_on_l1.bundleAt(lev);

    crescent::nested::bundles::interpolateParentBundleToChild(
        coupling_derivatives.core_snapshots.np1.bundleAt(0),
        level_context.coreLevelGeometry(),
        level_context.levelGeometry(),
        level_context.coreToL1RefinementRatio(),
        nested_config.derivative_interpolation_order,
        coupling_derivatives.core_interp_patch_scratch,
        parent_on_l1);
    crescent::nested::bundles::blendBundleWithParent(
        native_l1_bundle,
        parent_on_l1,
        level_context.nestedRuntime().coupling_zones.sigma_blend[lev],
        coupling_derivatives.l1_blended_derivatives.bundleAt(lev));
}

} // namespace crescent::nested::derivative_coupling

#endif // CRESCENT_NESTED_DERIVATIVE_COUPLING_H

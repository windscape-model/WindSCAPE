// Nested L1 post-step derivative recomputation for output fields.

#include "CRESCENT_NestedDerivativeCoupling.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedFreeSurfaceEvolution.H"
#include "CRESCENT_NestedPeriodicSurfaceDerivatives.H"
#include "CRESCENT_NestedRuntimeLogging.H"
#include "CRESCENT_NestedTypes.H"
#include "shared/CRESCENT_LevelEvolutionContext.H"

namespace {

using LevelEvolutionContext = crescent::shared::LevelEvolutionContext;
namespace derivative_coupling = crescent::nested::derivative_coupling;
namespace periodic_surface_derivatives =
    crescent::nested::periodic_surface_derivatives;

} // namespace

namespace crescent::nested::free_surface_evolution {

void recomputeL1OutputDerivatives(LevelEvolutionContext& level_context)
{
    const int lev = level_context.level();
    if (!crescent::nested::isActiveNestedL1Level(
            level_context.nestedConfig(),
            lev)) {
        return;
    }

    auto& coupling_derivatives =
        level_context.nestedRuntime().coupling_derivatives;
    const int verbosity = level_context.runtime().verbosity;
    crescent::nested::runtime_logging::logNestedOutputDerivativeRecomputeStart(
        verbosity);

    auto native_l1_bundle =
        level_context.solutionState().derivativeBundleAt(lev);
    const auto output_views =
        crescent::nested::bundles::makeL1OutputDerivativeViews(
            coupling_derivatives,
            lev);
    const auto final_output_blend_context =
        derivative_coupling::makeNestedFinalOutputBlendContext(
            level_context,
            native_l1_bundle);

    auto& periodic_fields =
        level_context.nestedRuntime().periodic_derivative_scratch[lev];
    periodic_surface_derivatives::prepareFinalOutputPeriodicL1Derivatives(
        level_context,
        periodic_fields);
    derivative_coupling::prepareNestedL1FinalOutputDerivatives(
        final_output_blend_context,
        output_views.core_np1,
        output_views.core_derivatives_on_l1,
        output_views.output_derivatives);

    crescent::nested::runtime_logging::logNestedOutputDerivativeRecomputeComplete(
        verbosity);
}

} // namespace crescent::nested::free_surface_evolution

#ifndef CRESCENT_NESTED_PERIODIC_SURFACE_DERIVATIVES_H
#define CRESCENT_NESTED_PERIODIC_SURFACE_DERIVATIVES_H

namespace crescent::nested {

struct PeriodicDerivativeScratch;

namespace stage_diagnostics {
struct RK4StageDiagnosticsContext;
}

} // namespace crescent::nested

namespace crescent::shared {
class LevelEvolutionContext;
}

namespace crescent::nested::periodic_surface_derivatives {

void prepareRK4StagePeriodicL1Derivatives(
    crescent::shared::LevelEvolutionContext& level_context,
    PeriodicDerivativeScratch& periodic_fields,
    stage_diagnostics::RK4StageDiagnosticsContext& stage_diagnostics);

void prepareFinalOutputPeriodicL1Derivatives(
    crescent::shared::LevelEvolutionContext& level_context,
    PeriodicDerivativeScratch& periodic_fields);

} // namespace crescent::nested::periodic_surface_derivatives

#endif // CRESCENT_NESTED_PERIODIC_SURFACE_DERIVATIVES_H

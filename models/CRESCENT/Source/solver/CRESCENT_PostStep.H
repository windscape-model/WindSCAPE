#ifndef CRESCENT_POST_STEP_H
#define CRESCENT_POST_STEP_H

namespace crescent::shared {

class LevelEvolutionContext;

} // namespace crescent::shared

namespace crescent::solver::poststep {

void applyOptionalPostTimestepLowPassFilter(
    const crescent::shared::LevelEvolutionContext& level_context);

void recomputeStandardDerivativeBundle(
    const crescent::shared::LevelEvolutionContext& level_context);

void computeVelocityOutputs(
    const crescent::shared::LevelEvolutionContext& level_context);

} // namespace crescent::solver::poststep

#endif // CRESCENT_POST_STEP_H

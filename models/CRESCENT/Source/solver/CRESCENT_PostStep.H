#ifndef CRESCENT_POST_STEP_H
#define CRESCENT_POST_STEP_H

namespace crescent::shared {
class LevelContext;
}

namespace crescent::solver::poststep {

void recomputeStandardDerivativeBundle(
    const crescent::shared::LevelContext& level_context);

void computeVelocityOutputs(
    const crescent::shared::LevelContext& level_context);

} // namespace crescent::solver::poststep

#endif // CRESCENT_POST_STEP_H

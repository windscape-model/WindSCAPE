#ifndef CRESCENT_NESTED_TOTAL_FIELD_H
#define CRESCENT_NESTED_TOTAL_FIELD_H

// Periodized total-field coupling (PTF, crescent.l1_coupling_mode = 1): the
// child tendency from C3-repaired periodic copies of the child state, blended
// with the interpolated parent tendency, and the post-step state relaxation.

namespace crescent::shared {
class LevelContext;
}

namespace crescent::nested::total_field {

// RHS of RK4 stage irk: (1 - sigma) F_L1(periodic copies) + sigma I(F_L0)
// at the same parent stage.
void runStage(
    crescent::shared::LevelContext& level_context,
    int irk);

// Post-RK4 state relaxation u <- (1 - sigma) P[u] + sigma I(u_L0), where
// P[u] is the periodic copy of the updated child state.
void applyStateRelaxation(
    crescent::shared::LevelContext& level_context);

// Output bundle of the current child state: periodic copies, their
// derivatives and vertical velocity, blended with the parent bundle. The
// preparation is kept for stage 0 of the next step.
void prepareOutputDerivatives(
    crescent::shared::LevelContext& level_context);

} // namespace crescent::nested::total_field

#endif // CRESCENT_NESTED_TOTAL_FIELD_H

#ifndef CRESCENT_NESTED_CORRECTION_FIELD_H
#define CRESCENT_NESTED_CORRECTION_FIELD_H

// Correction-field (R2) coupling, l1_coupling_mode == 2.
// Spec: docs/R2_CORRECTION_FIELD_SPEC.md (v1.1, adversarially verified).
// The child evolves TOTAL = interpolated-parent background + compact
// correction; every child FFT input is either compactly supported or an
// exact finite Fourier series. No C3 repair anywhere in this path.

#include <AMReX_REAL.H>

namespace crescent::shared {
class LevelContext;
}

namespace crescent::nested::correction_field {

// Observational capture at the end of each parent RK stage (copies + four
// parent-grid spectral quantities). The parent's own evolution is untouched.
void captureCoreStage(
    const crescent::shared::LevelContext& core_level_context,
    int rk_stage);

// One child RK stage: baseline reconstruction, tapered correction, delta
// derivatives/DNO + order-2 cross terms, DeltaF assembly, RHS write.
void runL1CorrectionFieldStage(
    crescent::shared::LevelContext& level_context,
    int irk,
    amrex::Real stage_time,
    amrex::Real temporal_weight);

// Post-RK4: u <- I(u_L0^{n+1}) + HL[chi * (u - I(u_L0^{n+1}))].
void applyPostRK4CorrectionFieldTaper(
    crescent::shared::LevelContext& level_context);

// Output derivative bundle for mode 2: interpolated parent np1 bundle plus
// spectral derivatives/DNO of the final (filtered, Fourier-series)
// correction, written into l1_blended_derivatives storage.
void recomputeCorrectionFieldOutputDerivatives(
    crescent::shared::LevelContext& level_context);

} // namespace crescent::nested::correction_field

#endif // CRESCENT_NESTED_CORRECTION_FIELD_H

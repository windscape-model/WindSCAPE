#ifndef CRESCENT_NESTED_RK4_STAGE_H
#define CRESCENT_NESTED_RK4_STAGE_H

#include "CRESCENT_NestedFreeSurfaceEvolution.H"

namespace crescent::nested::rk4_stage {

void runL1RK4Stage(
    crescent::shared::LevelEvolutionContext& level_context,
    int irk,
    amrex::Real stage_time,
    amrex::Real temporal_weight);

} // namespace crescent::nested::rk4_stage

#endif // CRESCENT_NESTED_RK4_STAGE_H

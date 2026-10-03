#pragma once

#include "CRESCENT_LevelInitializationContext.H"

namespace crescent::solver::setup::derived_fields {

void initializeCoreDerivedFields(
    const crescent::shared::LevelInitializationContext& level_context);

} // namespace crescent::solver::setup::derived_fields

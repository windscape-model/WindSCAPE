#pragma once

#include "CRESCENT_NestedConfig.H"

namespace amrex {

class ParmParse;

} // namespace amrex

namespace crescent::shared {

struct WaveRuntimeState;

} // namespace crescent::shared

namespace crescent::nested::config {

struct UserSpecifiedOptions
{
    bool blending_steepness = false;
    bool l1_refinement_box = false;
};

struct ParsedAlgorithmInputs
{
    AlgorithmConfig config;
    UserSpecifiedOptions user_specified_options;
};

ParsedAlgorithmInputs readAlgorithmInputs(
    amrex::ParmParse& pp);

void validateAlgorithmSelection(
    const AlgorithmConfig& nested_config,
    int max_level);

void validateRefinementRegionInputs(
    const AlgorithmConfig& nested_config,
    const UserSpecifiedOptions& user_specified_options,
    int max_level);

void validateInterpolationInputs(const AlgorithmConfig& nested_config);

void validateBlendingInputs(const AlgorithmConfig& nested_config);

void validateC3Inputs(const AlgorithmConfig& nested_config);

void emitInputWarnings(
    const AlgorithmConfig& nested_config,
    const UserSpecifiedOptions& user_specified_options);

} // namespace crescent::nested::config

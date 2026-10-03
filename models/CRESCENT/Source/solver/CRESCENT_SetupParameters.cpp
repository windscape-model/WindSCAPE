#include "CRESCENT_SetupParameters.H"

#include <utility>

#include <AMReX_ParmParse.H>

#include "CRESCENT_NestedAlgorithmInputs.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_ConstructionContext.H"
#include "CRESCENT_SetupRuntimeInputs.H"
#include "CRESCENT_WaveState.H"

namespace {

struct ConfiguredSetupState
{
    crescent::shared::RuntimeControlState runtime;
    crescent::shared::WaveRuntimeState wave;
    crescent::nested::AlgorithmConfig nested_config;
    crescent::nested::config::UserSpecifiedOptions nested_user_options;
};

ConfiguredSetupState buildConfiguredSetupState(
    int finest_level,
    const crescent::shared::RuntimeControlState& runtime,
    const crescent::shared::WaveRuntimeState& wave,
    const crescent::solver::setup::runtime_inputs::ParsedRuntimeInputs& runtime_inputs,
    const crescent::nested::config::ParsedAlgorithmInputs& nested_inputs)
{
    ConfiguredSetupState configured_state{
        runtime,
        wave,
        nested_inputs.config,
        nested_inputs.user_specified_options
    };
    crescent::solver::setup::runtime_inputs::applyRuntimeInputs(
        configured_state.runtime,
        configured_state.wave,
        finest_level,
        runtime_inputs);
    return configured_state;
}

void validateConfiguredSetupState(
    int max_level,
    const ConfiguredSetupState& configured_state)
{
    crescent::nested::config::validateAlgorithmSelection(
        configured_state.nested_config,
        max_level);
    crescent::nested::config::validateRefinementRegionInputs(
        configured_state.nested_config,
        configured_state.nested_user_options,
        max_level);
    crescent::nested::config::validateInterpolationInputs(
        configured_state.nested_config);
    crescent::solver::setup::runtime_inputs::validateWaveState(
        configured_state.wave);
    crescent::nested::config::validateBlendingInputs(
        configured_state.nested_config);
    crescent::nested::config::validateC3Inputs(
        configured_state.nested_config);
    crescent::nested::config::emitInputWarnings(
        configured_state.nested_config,
        configured_state.nested_user_options);
    crescent::solver::setup::runtime_inputs::validateRuntimeControlState(
        configured_state.runtime);
}

void applyConfiguredSetupState(
    crescent::shared::RuntimeControlState& runtime,
    crescent::shared::WaveRuntimeState& wave,
    crescent::nested::AlgorithmConfig& nested_config,
    ConfiguredSetupState configured_state)
{
    runtime = std::move(configured_state.runtime);
    wave = std::move(configured_state.wave);
    nested_config = std::move(configured_state.nested_config);
}

} // namespace

namespace crescent::solver::setup {

void readParameters(const ConstructionContext& construction_context)
{
    auto core = construction_context.coreAccess();
    auto nested = construction_context.nestedAccess();

    amrex::ParmParse crescent_pp(construction_context.solverParmParsePrefix());
    const auto runtime_inputs =
        crescent::solver::setup::runtime_inputs::readRuntimeInputs(
            crescent_pp);

    const auto nested_inputs =
        crescent::nested::config::readAlgorithmInputs(
            crescent_pp);
    auto configured_state = buildConfiguredSetupState(
        construction_context.finestLevel(),
        core.runtime,
        core.wave,
        runtime_inputs,
        nested_inputs);
    validateConfiguredSetupState(
        construction_context.maxLevel(),
        configured_state);
    applyConfiguredSetupState(
        core.runtime,
        core.wave,
        nested.config,
        std::move(configured_state));
}

} // namespace crescent::solver::setup

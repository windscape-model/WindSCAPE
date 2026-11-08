#pragma once

#include <iosfwd>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_WaveState.H"

namespace crescent::solver::checkpoint::config_io {

void writeWaveInitializationState(
    std::ostream& os,
    const crescent::shared::WaveInitializationState& initialization);

void readWaveInitializationState(
    std::istream& is,
    crescent::shared::WaveInitializationState& initialization);

void writeWaveSpectralState(
    std::ostream& os,
    const crescent::shared::WaveSpectralState& spectral);

void readWaveSpectralState(
    std::istream& is,
    crescent::shared::WaveSpectralState& spectral);

void writeNestedConfig(
    std::ostream& os,
    const crescent::nested::AlgorithmConfig& config);

void readNestedConfig(
    std::istream& is,
    crescent::nested::AlgorithmConfig& config);

} // namespace crescent::solver::checkpoint::config_io

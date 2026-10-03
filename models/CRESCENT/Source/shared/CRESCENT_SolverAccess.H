#ifndef CRESCENT_SOLVER_ACCESS_H
#define CRESCENT_SOLVER_ACCESS_H

#include <memory>

namespace crescent { class ProblemBase; }

namespace crescent::shared {

struct RootLayouts;
struct RuntimeControlState;
struct SolutionState;
struct SpectralWorkspace;
struct WaveRuntimeState;

} // namespace crescent::shared

namespace crescent::nested {

struct AlgorithmConfig;
struct RuntimeState;

} // namespace crescent::nested

namespace crescent::ownership {

struct CoreSolverAccess
{
    shared::RuntimeControlState& runtime;
    shared::WaveRuntimeState& wave;
    shared::SolutionState& solution;
    shared::RootLayouts& root_layouts;
    shared::SpectralWorkspace& spectral_workspace;
};

struct ConstCoreSolverAccess
{
    const shared::RuntimeControlState& runtime;
    const shared::WaveRuntimeState& wave;
    const shared::SolutionState& solution;
    const shared::RootLayouts& root_layouts;
    const shared::SpectralWorkspace& spectral_workspace;
};

struct NestedSolverAccess
{
    nested::AlgorithmConfig& config;
    nested::RuntimeState& runtime;
};

struct ConstNestedSolverAccess
{
    const nested::AlgorithmConfig& config;
    const nested::RuntimeState& runtime;
};

struct AmrIntegrationAccess
{
    std::unique_ptr<crescent::ProblemBase>& prob;
};

struct SolverAccess
{
    CoreSolverAccess core;
    NestedSolverAccess nested;
};

struct ConstSolverAccess
{
    ConstCoreSolverAccess core;
    ConstNestedSolverAccess nested;
};

} // namespace crescent::ownership

#endif // CRESCENT_SOLVER_ACCESS_H

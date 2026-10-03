#ifndef CRESCENT_SOLVER_OWNERSHIP_H
#define CRESCENT_SOLVER_OWNERSHIP_H

#include <memory>

#include "ProbBase.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedState.H"
#include "CRESCENT_RootLayouts.H"
#include "CRESCENT_RuntimeControl.H"
#include "CRESCENT_SolutionState.H"
#include "CRESCENT_SolverAccess.H"
#include "CRESCENT_SpectralWorkspace.H"
#include "CRESCENT_WaveState.H"

namespace crescent::ownership {

struct CoreSolverState
{
    shared::RuntimeControlState runtime;
    shared::WaveRuntimeState wave;
    shared::SolutionState solution;
    shared::RootLayouts root_layouts;
    shared::SpectralWorkspace spectral_workspace;
};

struct NestedSolverState
{
    nested::AlgorithmConfig config;
    nested::RuntimeState runtime;
};

struct AmrIntegrationState
{
    std::unique_ptr<crescent::ProblemBase> prob = nullptr;
};

struct SolverOwnership
{
    [[nodiscard]] CoreSolverAccess coreAccess() noexcept
    {
        return {
            core_.runtime,
            core_.wave,
            core_.solution,
            core_.root_layouts,
            core_.spectral_workspace
        };
    }

    [[nodiscard]] ConstCoreSolverAccess coreAccess() const noexcept
    {
        return {
            core_.runtime,
            core_.wave,
            core_.solution,
            core_.root_layouts,
            core_.spectral_workspace
        };
    }

    [[nodiscard]] NestedSolverAccess nestedAccess() noexcept
    {
        return {nested_.config, nested_.runtime};
    }

    [[nodiscard]] ConstNestedSolverAccess nestedAccess() const noexcept
    {
        return {nested_.config, nested_.runtime};
    }

    [[nodiscard]] SolverAccess solverAccess() noexcept
    {
        return {coreAccess(), nestedAccess()};
    }

    [[nodiscard]] ConstSolverAccess solverAccess() const noexcept
    {
        return {coreAccess(), nestedAccess()};
    }

    [[nodiscard]] AmrIntegrationAccess amrAccess() noexcept
    {
        return {amr_.prob};
    }

private:
    CoreSolverState core_;
    NestedSolverState nested_;
    AmrIntegrationState amr_;
};

} // namespace crescent::ownership

#endif // CRESCENT_SOLVER_OWNERSHIP_H

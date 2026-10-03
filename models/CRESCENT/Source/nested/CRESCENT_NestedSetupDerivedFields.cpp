#include "CRESCENT_NestedSetupDerivedFields.H"

#include <cmath>

#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_DNOExecution.H"
#include "CRESCENT_LevelWaveOperatorView.H"
#include "CRESCENT_Physics.H"

namespace {
void clearTinyInitializationOutput(amrex::MultiFab& mf)
{
    constexpr amrex::Real init_output_zero_tol = 1.0e-28;

    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        auto arr = mf.array(mfi);
        amrex::ParallelFor(
            bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (std::abs(arr(i, j, k)) < init_output_zero_tol) {
                    arr(i, j, k) = 0.0;
                }
            });
    }
}

} // namespace

namespace crescent::nested::setup::derived_fields {

void initializeFromState(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    const auto level_wave_operator =
        crescent::shared::makeLevelWaveOperatorView(
            lev,
            core,
            level_context.levelMesh().domain());
    auto derivative_bundle = core.solution.derivativeBundleAt(lev);
    auto& solution = core.solution;

    crescent::dno::computeSolutionStateDNO(
        level_wave_operator,
        core,
        solution.phi_z[lev]);

    if (core.runtime.verbosity >= 1) {
        amrex::Print() << "  L1 DNO computed: phi_z range = ["
                       << solution.phi_z[lev].min(0) << ", "
                       << solution.phi_z[lev].max(0) << "]" << std::endl;
    }

    crescent::physics::computeHorizontalVelocities(
        solution.u[lev],
        solution.v[lev],
        derivative_bundle,
        level_wave_operator.fft,
        level_wave_operator.geomdata,
        level_wave_operator.spectral_workspace.dealiasing_scratch[level_wave_operator.lev]);

    clearTinyInitializationOutput(solution.eta_x[lev]);
    clearTinyInitializationOutput(solution.eta_y[lev]);
    clearTinyInitializationOutput(solution.phi_x[lev]);
    clearTinyInitializationOutput(solution.phi_y[lev]);
    clearTinyInitializationOutput(solution.phi_z[lev]);
    clearTinyInitializationOutput(solution.u[lev]);
    clearTinyInitializationOutput(solution.v[lev]);

    if (core.runtime.verbosity >= 1) {
        amrex::Print() << "  L1 velocities computed: u range = ["
                       << solution.u[lev].min(0) << ", "
                       << solution.u[lev].max(0)
                       << "], v range = ["
                       << solution.v[lev].min(0) << ", "
                       << solution.v[lev].max(0) << "]" << std::endl;
    }
}

} // namespace crescent::nested::setup::derived_fields

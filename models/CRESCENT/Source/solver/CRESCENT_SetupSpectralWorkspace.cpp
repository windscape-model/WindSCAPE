#include "CRESCENT_SetupSpectralWorkspace.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

namespace {

void initializeLevelSpectralWorkspace(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    const auto core = level_context.coreAccess();
    auto& root = core.root_layouts;
    auto& workspace = core.spectral_workspace;
    const auto& wave = level_context.wave();
    const int verbosity = level_context.runtime().verbosity;
    const auto level_mesh = level_context.levelMesh();

    amrex::Box fft_box = level_mesh.grids.minimalBox();
    workspace.fft[lev] =
        std::make_unique<amrex::FFT::R2C<amrex::Real>>(fft_box);

    amrex::FFT::Info fft_info_root;
    fft_info_root.setNumProcs(1);
    workspace.fft_root[lev] =
        std::make_unique<amrex::FFT::R2C<amrex::Real>>(fft_box, fft_info_root);

    auto const& [cba, cdm] = workspace.fft[lev]->getSpectralDataLayout();
    if (verbosity >= 2) {
        amrex::Print() << " Domain BA lev " << lev << ": " << level_mesh.grids << std::endl;
        amrex::Print() << "   Spectral BA lev " << lev << ": " << cba << std::endl;
    }

    workspace.phi_expansion_spectrum[lev].define(
        cba, cdm, wave.spectral.perturbation_order, 0);

    workspace.compute_rhs_metric[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.compute_rhs_phi_z_sq[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.expansion_accum[lev].define(
        level_mesh.grids, level_mesh.dmap, wave.spectral.perturbation_order, 0);
    workspace.spectral_work[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.expansion_deriv[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.eta_spectral[lev].define(cba, cdm, 1, 0);
    workspace.deriv_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_input_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_x_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_y_spectral[lev].define(cba, cdm, 1, 0);
    workspace.directional_x_valid[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    workspace.directional_y_valid[lev].define(level_mesh.grids, level_mesh.dmap, 1, 0);
    // scalar_nghost follows the convention in SolutionState::defineScratchLevelStorage
    // so the scratch covers the widest nGrow that dealiasing consumers pass.
    const int dealiasing_nghost = (lev == 0) ? 1 : 2;
    workspace.dealiasing_scratch[lev].define(
        level_mesh.grids, level_mesh.dmap, 1, dealiasing_nghost);
    if (verbosity >= 2) {
        amrex::Print()
            << "   Allocated temp MultiFabs (real: metric, phi_z_sq, expansion_accum, spectral_work, expansion_deriv, directional_x_valid, directional_y_valid, dealiasing_scratch; spectral: eta_spectral, deriv_spectral, directional_input_spectral, directional_x_spectral, directional_y_spectral) for lev "
            << lev << std::endl;
    }

    core.solution.wavenumber_pow[lev].clear();
    core.solution.wavenumber_pow[lev].define(cba, cdm, wave.spectral.perturbation_order, 0);

    amrex::Vector<int> croot_map(1);
    croot_map[0] = amrex::ParallelDescriptor::IOProcessorNumber();
    root.cba[lev] = amrex::BoxArray(cba.minimalBox());
    root.cdmap[lev] = amrex::DistributionMapping(root.cba[lev]);
    root.cdmap[lev].define(croot_map);

    workspace.peak_wavenumber_gather_scratch[lev].define(
        root.cba[lev], root.cdmap[lev], 1, 0);
    workspace.peak_wavenumber_reduction_scratch[lev].define();
}

} // namespace

namespace crescent::solver::setup {

void initializeSpectralWorkspaces(
    const SolverInitializationContext& initialization_context)
{
    for (int lev = 0; lev <= initialization_context.finestLevel(); ++lev) {
        initializeLevelSpectralWorkspace(
            initialization_context.levelContext(lev));
    }
}

} // namespace crescent::solver::setup

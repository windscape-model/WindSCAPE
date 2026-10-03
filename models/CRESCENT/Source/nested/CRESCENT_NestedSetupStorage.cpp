#include "CRESCENT_NestedSetupStorage.H"

#include <AMReX_Print.H>

#include "CRESCENT_NestedC3Periodicity.H"
#include "CRESCENT_NestedDomains.H"

namespace {

void allocateNestedDerivativeCouplingStorage(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    auto nested = level_context.nestedAccess();
    const auto core_mesh = level_context.coreLevelMesh();
    const auto level_mesh = level_context.levelMesh();
    nested.runtime.coupling_derivatives.defineLevel(
        lev,
        core_mesh.grids, core_mesh.dmap,
        level_mesh.grids, level_mesh.dmap,
        level_context.refineRatio(),
        nested.config.blending_correction_y_smoothing_passes);
}

void allocateNestedPeriodicDerivativeScratch(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    auto core = level_context.coreAccess();
    auto nested = level_context.nestedAccess();
    const auto level_mesh = level_context.levelMesh();
    nested.runtime.periodic_derivative_scratch[lev].define(
        level_mesh.grids,
        level_mesh.dmap);

    // Post-RK4 state blend scratch. Matches the nghost of the level's
    // free-surface solution state (2 on L1+) so the scratch layout is
    // identical to the fresh-allocation pattern in applyPostRK4CoreStateBlend.
    const int post_rk4_nghost =
        core.solution.eta[lev].nGrow();
    nested.runtime
        .post_rk4_state_scratch[lev]
        .define(
            level_mesh.grids,
            level_mesh.dmap,
            post_rk4_nghost,
            level_context.refineRatio());

    // C3 repair gather workspace. The single-box global layout depends only on
    // the level's BoxArray.minimalBox(), which is fixed after initialization for
    // runs that do not regrid after setup.
    auto& c3_workspace =
        nested.runtime.c3_repair_workspace[lev];
    c3_workspace.define(level_mesh.grids.minimalBox());

    // Precompute the C3 interface basis. The basis depends only on
    // deterministic grid parameters, so building it once here avoids rebuilding
    // on every substage repair call.
    const amrex::Box global_box = level_mesh.grids.minimalBox();
    const int nx = global_box.length(0);
    const int ny = global_box.length(1);
    const int buffer_width =
        nested.runtime.coupling_zones.c3_buffer_width_computed[lev];
    const auto& level_geom = level_mesh.geom;
    crescent::nested::c3::computeC3InterfaceBasis(
        c3_workspace,
        crescent::nested::c3::InterfaceBasisOrder,
        nx,
        ny,
        buffer_width,
        buffer_width,
        level_geom.CellSize(0),
        level_geom.CellSize(1));
    c3_workspace.resizeLineScratches(
        crescent::nested::c3::InterfaceBasisOrder,
        nx,
        ny);

    nested.runtime.runtime_diagnostic_scratch[lev].define(
        level_mesh.grids,
        level_mesh.dmap,
        core.root_layouts.ba[lev],
        core.root_layouts.dmap[lev]);
}

} // namespace

namespace crescent::nested::setup::storage {

void initializeDomainsAndStorage(
    const crescent::shared::LevelInitializationContext& level_context)
{
    const int lev = level_context.level();
    auto core = level_context.coreAccess();
    auto nested = level_context.nestedAccess();
    const auto mesh = level_context.meshHierarchy();
    auto& runtime = core.runtime;
    const int verbosity = runtime.verbosity;
    if (verbosity >= 1) {
        amrex::Print() << "\n=== Setting up nested AMR buffer zone for L1 ===" << std::endl;
    }

    runtime.dt[lev] = runtime.dt[0];
    if (verbosity >= 1) {
        amrex::Print() << "Initialized dt[" << lev << "] = "
                       << runtime.dt[lev] << " s (same as L0)" << std::endl;
    }

    if (verbosity >= 1) {
        amrex::Print() << "\n=== Using CRESCENT Algorithm ===" << std::endl;
    }

    const auto coupling_context =
        crescent::nested::domains::makeCouplingZoneSetupContext(
            nested.config,
            nested.runtime,
            verbosity,
            mesh);
    crescent::nested::domains::initializeLevelCouplingZones(
        lev, coupling_context);
    allocateNestedDerivativeCouplingStorage(level_context);
    allocateNestedPeriodicDerivativeScratch(level_context);

    if (verbosity >= 1) {
        amrex::Print() << "=== CRESCENT initialization complete ===" << std::endl;
        amrex::Print() << "  Using L1-based FFT (fft[" << lev
                       << "]) for spectral operations" << std::endl;
    }
}

} // namespace crescent::nested::setup::storage

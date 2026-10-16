#include "CRESCENT_NestedSetup.H"

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include "CRESCENT_Interpolation.H"
#include "CRESCENT_NestedC3Periodicity.H"
#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedDomains.H"

namespace crescent::nested::setup {

namespace {

void allocateNestedDerivativeCouplingStorage(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& nested = level_context.nested();
    const auto core_mesh = level_context.coreLevelMesh();
    const auto level_mesh = level_context.levelMesh();
    nested.runtime.coupling_derivatives.defineLevel(
        lev,
        core_mesh.grids, core_mesh.dmap,
        level_mesh.grids, level_mesh.dmap,
        level_context.refineRatio());

    if (nested.config.l1_coupling_mode == 2) {
        nested.runtime.correction_field.defineCoreStages(
            core_mesh.grids, core_mesh.dmap);
        nested.runtime.correction_field.defineChild(
            level_mesh.grids, level_mesh.dmap);
    }
}

void allocateNestedPeriodicDerivativeScratch(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    auto& nested = level_context.nested();
    const auto level_mesh = level_context.levelMesh();
    nested.runtime.periodic_derivative_scratch[lev].define(
        level_mesh.grids,
        level_mesh.dmap);

    // Scratch of the post-RK4 state relaxation, with the nghost of the
    // level's free-surface state (2 on L1+).
    const int post_rk4_nghost =
        core.solution.eta[lev].nGrow();
    nested.runtime
        .post_rk4_state_scratch[lev]
        .define(
            level_mesh.grids,
            level_mesh.dmap,
            post_rk4_nghost,
            level_context.refineRatio());

    if (nested.config.l1_coupling_mode != 2) {
        // C3 repair gather workspace. CCF bypasses C3 completely, so avoid
        // constructing its global gather and basis caches in mode 2.
        auto& c3_workspace = nested.runtime.c3_repair_workspace[lev];
        c3_workspace.define(level_mesh.grids.minimalBox());

        const amrex::Box global_box = level_mesh.grids.minimalBox();
        const int nx = global_box.length(0);
        const int ny = global_box.length(1);
        const int buffer_width_x =
            nested.runtime.coupling_zones.c3_buffer_width_x_computed[lev];
        const int buffer_width_y =
            nested.runtime.coupling_zones.c3_buffer_width_y_computed[lev];
        const auto& level_geom = level_mesh.geom;
        crescent::nested::c3::computeC3InterfaceBasis(
            c3_workspace,
            crescent::nested::c3::InterfaceBasisOrder,
            nx,
            ny,
            buffer_width_x,
            buffer_width_y,
            level_geom.CellSize(0),
            level_geom.CellSize(1));
        c3_workspace.resizeLineScratches(nx, ny);
    }

    nested.runtime.runtime_diagnostic_scratch[lev].define(
        core.root_layouts.ba[lev],
        core.root_layouts.dmap[lev]);
}

void interpolateStateFromParent(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    AMREX_ALWAYS_ASSERT(lev > 0);
    auto& core = level_context.core();
    auto& nested = level_context.nested();
    const int verbosity = level_context.runtime().verbosity;
    if (verbosity >= 1) {
        amrex::Print() << "\n=== Initializing L1 from parent L0 ===" << std::endl;
    }

    const auto core_mesh = level_context.coreLevelMesh();
    const auto level_mesh = level_context.levelMesh();
    auto& solution = core.solution;
    auto& coarse_interp_patch =
        nested.runtime.post_rk4_state_scratch[lev].coarseInterpPatch();

    solution.eta[0].FillBoundary(core_mesh.geom.periodicity());
    solution.phi[0].FillBoundary(core_mesh.geom.periodicity());

    crescent::interpolation::interpolateParentToChild(
        solution.eta[0],
        solution.phi[0],
        solution.eta[lev],
        solution.phi[lev],
        core_mesh.geom,
        level_mesh.geom,
        level_context.refineRatio(),
        coarse_interp_patch,
        nested.config.state_interpolation_order);

    core.solution.saveSurfaceState(lev);
}

} // namespace

void initializeDomainsAndStorage(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    auto& nested = level_context.nested();
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
        crescent::nested::domains::CouplingZoneSetupContext{
            nested.config,
            nested.runtime,
            verbosity,
            mesh};
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

void zeroDerivedFields(
    const crescent::shared::LevelContext& level_context)
{
    // The live child state is a non-periodic window of the parent. Its
    // output bundle is assembled later from repaired periodic copies (PTF)
    // or the parent-plus-compact-correction split (CCF); never FFT it here.
    const int lev = level_context.level();
    auto& solution = level_context.core().solution;
    solution.eta_x[lev].setVal(0.0);
    solution.eta_y[lev].setVal(0.0);
    solution.phi_x[lev].setVal(0.0);
    solution.phi_y[lev].setVal(0.0);
    solution.phi_z[lev].setVal(0.0);
    solution.u[lev].setVal(0.0);
    solution.v[lev].setVal(0.0);
}

void initializeNestedL1Level(
    const crescent::shared::LevelContext& level_context)
{
    const int lev = level_context.level();
    auto& core = level_context.core();
    auto& nested = level_context.nested();
    AMREX_ALWAYS_ASSERT(crescent::nested::isActiveNestedL1Level(
        nested.config,
        lev));

    initializeDomainsAndStorage(level_context);
    interpolateStateFromParent(level_context);
    zeroDerivedFields(level_context);

    if (core.runtime.verbosity >= 1) {
        amrex::Print() << "=== Nested AMR buffer zone setup complete ===" << std::endl;
    }
}

} // namespace crescent::nested::setup

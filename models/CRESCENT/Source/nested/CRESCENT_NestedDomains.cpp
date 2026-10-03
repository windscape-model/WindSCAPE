#include "CRESCENT_NestedDomains.H"

#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>

#include "CRESCENT_NestedConfig.H"
#include "CRESCENT_NestedMasks.H"
#include "CRESCENT_NestedState.H"

namespace {

using CouplingZoneSetupContext =
    crescent::nested::domains::CouplingZoneSetupContext;
using L1RefinementRegion = crescent::nested::domains::L1RefinementRegion;
using L1RefinementTaggingContext =
    crescent::nested::domains::L1RefinementTaggingContext;

void initializeCouplingBufferWidth(int lev, const CouplingZoneSetupContext& context)
{
    const auto& nested = context.config;
    auto& coupling_zones = context.runtime.coupling_zones;
    const auto level_mesh = context.levelMesh(lev);
    const amrex::Real dx_nested = level_mesh.geom.CellSize(0);
    const amrex::Real distance_to_omega_S =
        coupling_zones.domains[lev].x_lo_S - coupling_zones.domains[lev].x_lo_L1;
    const int computed_buffer_width =
        static_cast<int>(std::ceil(distance_to_omega_S / dx_nested));

    coupling_zones.c3_buffer_width_computed[lev] = (nested.c3_use_interface == 1)
        ? std::max(computed_buffer_width, 8)
        : computed_buffer_width;

    if (context.verbosity >= 1) {
        amrex::Print() << "\n  C3-repair buffer width computation:" << std::endl;
        amrex::Print() << "    Distance to Omega_S: " << distance_to_omega_S << " m" << std::endl;
        amrex::Print() << "    Computed buffer: " << computed_buffer_width << " cells" << std::endl;
        amrex::Print() << "    Actual buffer:   "
                       << coupling_zones.c3_buffer_width_computed[lev] << " cells";
        if (nested.c3_use_interface == 1 &&
            coupling_zones.c3_buffer_width_computed[lev] > computed_buffer_width) {
            amrex::Print() << " (enforced minimum for 5-pt stencil)";
        }
        amrex::Print() << std::endl;
        if (coupling_zones.c3_buffer_width_computed[lev] > computed_buffer_width) {
            const amrex::Real extension =
                (coupling_zones.c3_buffer_width_computed[lev] - computed_buffer_width) * dx_nested;
            const amrex::Real x_end =
                coupling_zones.domains[lev].x_lo_L1 +
                coupling_zones.c3_buffer_width_computed[lev] * dx_nested;
            amrex::Print() << "    WARNING: C3-repair extends " << extension
                           << " m beyond Omega_S edge" << std::endl;
            amrex::Print() << "    (ends at x=" << x_end << " m, Omega_S starts at "
                           << coupling_zones.domains[lev].x_lo_S << " m)" << std::endl;
        }
    }
}

void initializeBlendingField(int lev, const CouplingZoneSetupContext& context)
{
    const auto& nested = context.config;
    auto& coupling_zones = context.runtime.coupling_zones;
    const auto level_mesh = context.levelMesh(lev);
    const int active_blending_mode =
        crescent::nested::resolveActiveBlendingMode(nested);

    if (active_blending_mode == 4) {
        if (context.verbosity >= 1) {
            amrex::Print() << "  Using harmonic 2D blending with septic endpoint smoothing" << std::endl;
        }
        crescent::nested::masks::generateBlendingFunctionHarmonicSmooth(
            coupling_zones.sigma_blend[lev],
            level_mesh.geom,
            coupling_zones.domains[lev],
            nested.blending_steepness);
    } else if (active_blending_mode == 5) {
        if (context.verbosity >= 1) {
            amrex::Print() << "  Using superellipse 2D blending, gamma="
                           << nested.blending_superellipse_gamma
                           << ", p=" << nested.blending_superellipse_p << std::endl;
        }
        crescent::nested::masks::generateBlendingFunctionSuperellipse(
            coupling_zones.sigma_blend[lev],
            level_mesh.geom,
            coupling_zones.domains[lev],
            nested.blending_superellipse_gamma,
            nested.blending_superellipse_p);
    } else if (active_blending_mode == 3) {
        if (context.verbosity >= 1) {
            amrex::Print() << "  Using tensor-product 2D blending (smooth 1D ramps with corner composition), alpha="
                           << nested.blending_tensor_alpha << std::endl;
        }
        crescent::nested::masks::generateBlendingFunctionTensorProduct(
            coupling_zones.sigma_blend[lev],
            level_mesh.geom,
            coupling_zones.domains[lev],
            nested.blending_tensor_alpha);
    } else if (active_blending_mode == 2) {
        if (context.verbosity >= 1) {
            amrex::Print() << "  Using harmonic 2D blending (Omega_A=0, Omega_S boundary=1)" << std::endl;
        }
        crescent::nested::masks::generateBlendingFunctionHarmonic(
            coupling_zones.sigma_blend[lev],
            level_mesh.geom,
            coupling_zones.domains[lev],
            nested.blending_steepness);
    } else if (active_blending_mode == 1) {
        if (context.verbosity >= 1) {
            amrex::Print() << "  Using sequential 1D blending (preserves y-uniformity)" << std::endl;
        }
        crescent::nested::masks::generateBlendingFunctionSequential(
            coupling_zones.sigma_blend[lev],
            level_mesh.geom,
            coupling_zones.domains[lev],
            nested.blending_steepness);
    } else {
        if (context.verbosity >= 1) {
            amrex::Print() << "  Using 2D Euclidean blending" << std::endl;
        }
        crescent::nested::masks::generateBlendingFunction(
            coupling_zones.sigma_blend[lev],
            level_mesh.geom,
            coupling_zones.domains[lev],
            nested.blending_steepness);
    }
}

bool onBaseLevelEdge(amrex::Real val, amrex::Real dx, amrex::Real origin)
{
    const amrex::Real offset = val - origin;
    return std::abs(std::round(offset / dx) * dx - offset) < 1.0e-10;
}

L1RefinementRegion readL1RefinementRegion(
    const L1RefinementTaggingContext& context)
{
    const amrex::Real* prob_lo = context.geom.ProbLo();
    const amrex::Real* dx_L0 = context.geom.CellSize();

    L1RefinementRegion region{
        context.config.l1_box_x_lo,
        context.config.l1_box_x_hi,
        context.config.l1_box_y_lo,
        context.config.l1_box_y_hi,
        dx_L0[0] / context.ref_ratio[0],
        dx_L0[1] / context.ref_ratio[1],
        context.config.l1_box_x_hi - context.config.l1_box_x_lo,
        context.config.l1_box_y_hi - context.config.l1_box_y_lo,
        static_cast<int>(std::round(
            (context.config.l1_box_x_hi - context.config.l1_box_x_lo) /
            (dx_L0[0] / context.ref_ratio[0]))),
        static_cast<int>(std::round(
            (context.config.l1_box_y_hi - context.config.l1_box_y_lo) /
            (dx_L0[1] / context.ref_ratio[1])))
    };

    if (!onBaseLevelEdge(region.x_lo, dx_L0[0], prob_lo[0]) ||
        !onBaseLevelEdge(region.x_hi, dx_L0[0], prob_lo[0]) ||
        !onBaseLevelEdge(region.y_lo, dx_L0[1], prob_lo[1]) ||
        !onBaseLevelEdge(region.y_hi, dx_L0[1], prob_lo[1])) {
        amrex::Print() << "\n*** ERROR: L1 boundaries must fall on L0 cell edges ***\n"
                       << "  crescent.L1_box_lo = [" << region.x_lo << ", " << region.y_lo << "]\n"
                       << "  crescent.L1_box_hi = [" << region.x_hi << ", " << region.y_hi << "]\n"
                       << "  dx_L0 = " << dx_L0[0] << " m\n"
                       << "  Boundaries must be aligned with L0 cell edges measured from geometry.prob_lo.\n";
        amrex::Abort("L1 domain boundaries not on L0 cell edges.");
    }

    if (region.nx % context.blocking_factor[0] != 0 ||
        region.ny % context.blocking_factor[1] != 0) {
        amrex::Print() << "\n*** ERROR: L1 cell count must be a multiple of blocking_factor ***\n"
                       << "  L1 cells: " << region.nx << " x " << region.ny << "\n"
                       << "  blocking_factor: " << context.blocking_factor[0] << " x "
                       << context.blocking_factor[1] << "\n"
                       << "  L1 extent must be a multiple of "
                       << context.blocking_factor[0] * region.dx_x
                       << " m (blocking_factor * dx_L1)\n";
        amrex::Abort("L1 cell count not a multiple of blocking_factor.");
    }

    return region;
}

void printL1RefinementSummary(
    const L1RefinementTaggingContext& context,
    const L1RefinementRegion& region)
{
    const amrex::Real* prob_lo = context.geom.ProbLo();
    const amrex::Real* prob_hi = context.geom.ProbHi();
    const amrex::Real* dx_L0 = context.geom.CellSize();

    amrex::Print() << "\n=== Tagging cells for refinement (L0 -> L1) ===" << std::endl;
    amrex::Print() << "  L0 domain: [" << prob_lo[0] << ", " << prob_hi[0]
                   << "] x [" << prob_lo[1] << ", " << prob_hi[1]
                   << "], dx=" << dx_L0[0] << " m" << std::endl;
    amrex::Print() << "  L1 domain: [" << region.x_lo << ", " << region.x_hi
                   << "] x [" << region.y_lo << ", " << region.y_hi
                   << "], dx=" << region.dx_x << " m" << std::endl;
    amrex::Print() << "  L1 grid:   " << region.nx << " x " << region.ny
                   << " cells (" << region.Lx << " x " << region.Ly << " m)" << std::endl;
}

void applyL1RefinementTags(
    amrex::TagBoxArray& tags,
    const L1RefinementTaggingContext& context,
    const L1RefinementRegion& region)
{
    const amrex::Real* prob_lo = context.geom.ProbLo();
    const amrex::Real* dx_L0 = context.geom.CellSize();

    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        auto tag_arr = tags.array(mfi);

        const amrex::Real d_prob_lo_0 = prob_lo[0];
        const amrex::Real d_prob_lo_1 = prob_lo[1];
        const amrex::Real d_dx_0 = dx_L0[0];
        const amrex::Real d_dx_1 = dx_L0[1];
        const amrex::Real d_L1_lo_0 = region.x_lo;
        const amrex::Real d_L1_hi_0 = region.x_hi;
        const amrex::Real d_L1_lo_1 = region.y_lo;
        const amrex::Real d_L1_hi_1 = region.y_hi;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const amrex::Real x = d_prob_lo_0 + (i + 0.5) * d_dx_0;
            const amrex::Real y = d_prob_lo_1 + (j + 0.5) * d_dx_1;

            if (x >= d_L1_lo_0 && x < d_L1_hi_0 &&
                y >= d_L1_lo_1 && y < d_L1_hi_1) {
                tag_arr(i,j,k) = amrex::TagBox::SET;
            }
        });
    }
}

} // namespace

namespace crescent::nested::domains {

void initializeLevelCouplingZones(
    int lev,
    const CouplingZoneSetupContext& context)
{
    auto& coupling_zones = context.runtime.coupling_zones;
    const auto level_mesh = context.levelMesh(lev);
    coupling_zones.defineLevel(lev, level_mesh.grids, level_mesh.dmap);

    coupling_zones.domains[lev] = computeDomainStructure(
        level_mesh.geom,
        level_mesh.grids,
        context.config.roi_accurate_fraction,
        context.config.roi_spectral_fraction);

    if (context.verbosity >= 1) {
        printDomainStructure(coupling_zones.domains[lev]);
    }
    initializeCouplingBufferWidth(lev, context);
    initializeBlendingField(lev, context);
}

void tagL1RefinementRegion(
    amrex::TagBoxArray& tags,
    const L1RefinementTaggingContext& context)
{
    const auto region = readL1RefinementRegion(context);
    if (context.verbosity >= 1) {
        printL1RefinementSummary(context, region);
    }
    applyL1RefinementTags(tags, context, region);
    if (context.verbosity >= 1) {
        amrex::Print() << "=== Tagging complete ===" << std::endl;
    }
}

DomainStructure computeDomainStructure(
    const amrex::Geometry& geom,
    const amrex::BoxArray& grids,
    amrex::Real accurate_fraction,
    amrex::Real spectral_fraction)
{
    using namespace amrex;

    DomainStructure domain;

    // Store fractions
    domain.accurate_fraction = accurate_fraction;
    domain.spectral_fraction = spectral_fraction;

    // Get actual L1 refined box extent (minimal bounding box of all grids)
    const Box& grid_box = grids.minimalBox();
    const Real* dx = geom.CellSize();
    const auto& prob_lo_global = geom.ProbLoArray();

    // Convert grid_box to physical coordinates
    // CRITICAL: This gives us the ACTUAL L1 refined region, not the logical domain
    Real actual_lo_x = prob_lo_global[0] + grid_box.smallEnd(0) * dx[0];
    Real actual_hi_x = prob_lo_global[0] + (grid_box.bigEnd(0) + 1) * dx[0];
    Real actual_lo_y = prob_lo_global[1] + grid_box.smallEnd(1) * dx[1];
    Real actual_hi_y = prob_lo_global[1] + (grid_box.bigEnd(1) + 1) * dx[1];

    // Store L1 box boundaries
    domain.x_lo_L1 = actual_lo_x;
    domain.x_hi_L1 = actual_hi_x;
    domain.y_lo_L1 = actual_lo_y;
    domain.y_hi_L1 = actual_hi_y;

    // L1 actual extent
    domain.Lx = actual_hi_x - actual_lo_x;
    domain.Ly = actual_hi_y - actual_lo_y;
    domain.ni = grid_box.length(0);
    domain.nj = grid_box.length(1);

    // Compute margins RELATIVE TO ACTUAL L1 EXTENT
    Real margin_A_x = accurate_fraction * domain.Lx;
    Real margin_A_y = accurate_fraction * domain.Ly;
    Real margin_S_x = spectral_fraction * domain.Lx;
    Real margin_S_y = spectral_fraction * domain.Ly;

    // Accurate domain Ω_A: centered within L1 box with specified margin
    domain.x_lo_A = actual_lo_x + margin_A_x;
    domain.x_hi_A = actual_hi_x - margin_A_x;
    domain.y_lo_A = actual_lo_y + margin_A_y;
    domain.y_hi_A = actual_hi_y - margin_A_y;

    // Spectral subdomain Ω_S: centered within L1 box with specified margin
    domain.x_lo_S = actual_lo_x + margin_S_x;
    domain.x_hi_S = actual_hi_x - margin_S_x;
    domain.y_lo_S = actual_lo_y + margin_S_y;
    domain.y_hi_S = actual_hi_y - margin_S_y;

    // Sanity check: Ω_A must be inside Ω_S
    AMREX_ALWAYS_ASSERT(domain.x_lo_A >= domain.x_lo_S);
    AMREX_ALWAYS_ASSERT(domain.x_hi_A <= domain.x_hi_S);
    AMREX_ALWAYS_ASSERT(domain.y_lo_A >= domain.y_lo_S);
    AMREX_ALWAYS_ASSERT(domain.y_hi_A <= domain.y_hi_S);

    // Compute zone widths
    domain.W_exterior_x = margin_S_x;  // Width of exterior zone Ω_E [x_lo_L1, x_lo_S]
    domain.W_exterior_y = margin_S_y;
    domain.W_mix_x = margin_A_x - margin_S_x;  // Width of mixing zone Ω_M [x_lo_S, x_lo_A]
    domain.W_mix_y = margin_A_y - margin_S_y;

    AMREX_ALWAYS_ASSERT(domain.W_mix_x > 0);
    AMREX_ALWAYS_ASSERT(domain.W_mix_y > 0);

    // Convert to index space
    // CRITICAL: Compute indices relative to GLOBAL coordinate system

    // Compute cell indices (round to nearest, relative to global origin)
    int i_lo_A = static_cast<int>(std::round((domain.x_lo_A - prob_lo_global[0]) / dx[0]));
    int i_hi_A = static_cast<int>(std::round((domain.x_hi_A - prob_lo_global[0]) / dx[0])) - 1;
    int j_lo_A = static_cast<int>(std::round((domain.y_lo_A - prob_lo_global[1]) / dx[1]));
    int j_hi_A = static_cast<int>(std::round((domain.y_hi_A - prob_lo_global[1]) / dx[1])) - 1;

    int i_lo_S = static_cast<int>(std::round((domain.x_lo_S - prob_lo_global[0]) / dx[0]));
    int i_hi_S = static_cast<int>(std::round((domain.x_hi_S - prob_lo_global[0]) / dx[0])) - 1;
    int j_lo_S = static_cast<int>(std::round((domain.y_lo_S - prob_lo_global[1]) / dx[1]));
    int j_hi_S = static_cast<int>(std::round((domain.y_hi_S - prob_lo_global[1]) / dx[1])) - 1;

    // Create boxes (k dimension is 0:0 for 2D surface)
    domain.box_A = Box(IntVect(i_lo_A, j_lo_A, 0),
                       IntVect(i_hi_A, j_hi_A, 0));
    domain.box_S = Box(IntVect(i_lo_S, j_lo_S, 0),
                       IntVect(i_hi_S, j_hi_S, 0));

    return domain;
}

void printDomainStructure(const DomainStructure& domain)
{
    using namespace amrex;

    Print() << "\n";
    Print() << "========================================\n";
    Print() << "CRESCENT Domain Structure\n";
    Print() << "========================================\n";
    Print() << "\n";
    Print() << "Configuration:\n";
    Print() << "  Accurate margin (Ω_A): " << domain.accurate_fraction * 100 << "%\n";
    Print() << "  Spectral margin (Ω_S): " << domain.spectral_fraction * 100 << "%\n";
    Print() << "  Grid: " << domain.ni << " x " << domain.nj << "\n";
    Print() << "  Domain size: Lx=" << domain.Lx << " m, Ly=" << domain.Ly << " m\n";
    Print() << "\n";

    Print() << "Physical Boundaries:\n";
    Print() << "  Ω_A (Accurate):  x=[" << domain.x_lo_A << ", " << domain.x_hi_A << "] m\n";
    Print() << "                   y=[" << domain.y_lo_A << ", " << domain.y_hi_A << "] m\n";
    Print() << "  Ω_S (Spectral):  x=[" << domain.x_lo_S << ", " << domain.x_hi_S << "] m\n";
    Print() << "                   y=[" << domain.y_lo_S << ", " << domain.y_hi_S << "] m\n";
    Print() << "\n";

    Print() << "Zone Widths:\n";
    Print() << "  Exterior zone (Ω_E):  W_exterior_x=" << domain.W_exterior_x << " m, W_exterior_y=" << domain.W_exterior_y << " m\n";
    Print() << "  Mixing zone (Ω_M):    W_mix_x=" << domain.W_mix_x << " m, W_mix_y=" << domain.W_mix_y << " m\n";
    Print() << "\n";

    Print() << "Index Space Boundaries:\n";
    Print() << "  Ω_A box: " << domain.box_A << "\n";
    Print() << "  Ω_S box: " << domain.box_S << "\n";
    Print() << "\n";

    // Compute coverage percentages
    Real coverage_A_x = (domain.x_hi_A - domain.x_lo_A) / domain.Lx * 100;
    Real coverage_A_y = (domain.y_hi_A - domain.y_lo_A) / domain.Ly * 100;
    Real coverage_S_x = (domain.x_hi_S - domain.x_lo_S) / domain.Lx * 100;
    Real coverage_S_y = (domain.y_hi_S - domain.y_lo_S) / domain.Ly * 100;

    Print() << "Coverage:\n";
    Print() << "  Ω_A: " << coverage_A_x << "% x " << coverage_A_y << "% of domain\n";
    Print() << "  Ω_S: " << coverage_S_x << "% x " << coverage_S_y << "% of domain\n";
    Print() << "  Ω_M: Ω_S \\ Ω_A (mixing zone)\n";
    Print() << "  Ω_E: L1 \\ Ω_S (exterior zone)\n";
    Print() << "\n";
    Print() << "========================================\n";
    Print() << "\n";
}

} // namespace crescent::nested::domains

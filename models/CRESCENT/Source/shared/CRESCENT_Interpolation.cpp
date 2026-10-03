#include "CRESCENT_Interpolation.H"

#include <AMReX_GpuLaunch.H>

namespace crescent::interpolation {

using amrex::Box;
using amrex::Geometry;
using amrex::IndexType;
using amrex::IntVect;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::Real;
using amrex::TilingIfNotGPU;
using amrex::Vector;

// ============================================================================
// Custom 2D quadratic interpolation for nz=1 (3D build, 2D physics)
// ============================================================================

// 1D quadratic Lagrange weights on coarse-cell centers at x = {-1,0,+1}
// evaluated at s in coarse-cell units (for ratio=2, s = +/- 0.25).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void quadraticLagrangeWeights(Real s, Real& wm1, Real& w0, Real& wp1) noexcept
{
    // L_{-1}(s), L_0(s), L_{+1}(s) for nodes -1,0,+1
    wm1 = Real(0.5) * s * (s - Real(1.0));
    w0  = Real(1.0) - s * s;
    wp1 = Real(0.5) * s * (s + Real(1.0));
}

void fillCoarseInterpolationPatch(
    const MultiFab& parent,
    MultiFab& coarse_patch,
    const Geometry& parent_geom,
    int ng_crse)
{
    AMREX_ALWAYS_ASSERT(coarse_patch.nGrow() >= ng_crse);

    coarse_patch.setVal(0.0);
    coarse_patch.ParallelCopy(
        parent,
        0,
        0,
        1,
        IntVect(0),
        IntVect(ng_crse, ng_crse, 0),
        parent_geom.periodicity());
    coarse_patch.FillBoundary(parent_geom.periodicity());
}

// Prolong from coarse to fine using 2D quadratic interpolation in x/y only.
// - Assumes ref_ratio = (2,2,1) and nz=1 (k=0 plane).
// - Fills child valid + ng_fine ghosts.
// - Requires parent plane (k=0) valid (ghosts optional; we build our own patch).
void prolongQuadraticRef2XY(
    const MultiFab& parent,        // coarse (L0), cell-centered, k=0 plane
    MultiFab&       child,         // fine (L1), cell-centered, k=0 plane, has ng_fine ghosts
    const Geometry& parent_geom,
    const Geometry& child_geom,
    const IntVect&  ref_ratio,
    int             ng_fine,        // fine ghosts to fill (must be <= child.nGrow())
    int             ng_crse,        // coarse patch ghosts (>=1 for quadratic; recommend 2)
    MultiFab&       crse_patch)
{
    AMREX_ALWAYS_ASSERT(ref_ratio[0] == 2 && ref_ratio[1] == 2 && ref_ratio[2] == 1);
    AMREX_ALWAYS_ASSERT(parent.ixType() == IndexType::TheCellType());
    AMREX_ALWAYS_ASSERT(child.ixType()  == IndexType::TheCellType());
    AMREX_ALWAYS_ASSERT(child.nGrow() >= ng_fine);

    // Keep the coarse valid region non-overlapping and let ParallelCopy fill
    // the coarse ghost cells directly from the parent. This preserves the
    // original interpolation stencil support without creating a rank-0-only
    // global patch.
    fillCoarseInterpolationPatch(parent, crse_patch, parent_geom, ng_crse);

    // ----- Step 2: Fill fine child (valid+ghost) using our own x/y quadratic kernel -----

    for (MFIter mfi(child, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Region to fill: valid + ng_fine ghost (clipped to domain)
        Box bx = mfi.validbox();
        bx.grow(IntVect(ng_fine, ng_fine, 0));
        bx &= child_geom.Domain();

        auto const& a_fine = child.array(mfi);
        auto const& a_crse = crse_patch.const_array(mfi);

        // Quadratic weights for the two fine positions inside coarse cell at ratio=2:
        // fine i even -> s = -0.25, fine i odd -> s = +0.25
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // nz=1 assumed
            if (k != 0) return;

            // coarse indices (ratio 2)
            const int ic = i >> 1; // i/2
            const int jc = j >> 1;

            // subcell offset s in coarse units
            const Real sx = ((i & 1) == 0) ? Real(-0.25) : Real(0.25);
            const Real sy = ((j & 1) == 0) ? Real(-0.25) : Real(0.25);

            Real wxm1, wx0, wxp1;
            Real wym1, wy0, wyp1;
            quadraticLagrangeWeights(sx, wxm1, wx0, wxp1);
            quadraticLagrangeWeights(sy, wym1, wy0, wyp1);

            // 3x3 tensor product stencil in x/y; k fixed = 0
            Real v = Real(0.0);

            // dy = -1
            v += wxm1 * wym1 * a_crse(ic-1, jc-1, 0);
            v += wx0  * wym1 * a_crse(ic  , jc-1, 0);
            v += wxp1 * wym1 * a_crse(ic+1, jc-1, 0);

            // dy = 0
            v += wxm1 * wy0  * a_crse(ic-1, jc  , 0);
            v += wx0  * wy0  * a_crse(ic  , jc  , 0);
            v += wxp1 * wy0  * a_crse(ic+1, jc  , 0);

            // dy = +1
            v += wxm1 * wyp1 * a_crse(ic-1, jc+1, 0);
            v += wx0  * wyp1 * a_crse(ic  , jc+1, 0);
            v += wxp1 * wyp1 * a_crse(ic+1, jc+1, 0);

            a_fine(i,j,0) = v;
        });
    }

    // ----- Step 3: Fill fine ghosts across boxes/ranks (periodic) -----
    child.FillBoundary(child_geom.periodicity());
}

// Bilinear interpolation (1st order, simpler, works without fine ghosts)
void prolongBilinearRef2XY(
    const MultiFab& parent,
    MultiFab&       child,
    const Geometry& parent_geom,
    const Geometry& child_geom,
    const IntVect&  ref_ratio,
    MultiFab&       crse_patch)
{
    AMREX_ALWAYS_ASSERT(ref_ratio[0] == 2 && ref_ratio[1] == 2 && ref_ratio[2] == 1);

    // Bilinear only needs valid region access (no fine ghosts required)
    const int ng_crse = 1; // minimal coarse ghost for periodic BC

    fillCoarseInterpolationPatch(parent, crse_patch, parent_geom, ng_crse);

    // Bilinear interpolation kernel
    for (MFIter mfi(child, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto const& a_fine = child.array(mfi);
        auto const& a_crse = crse_patch.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k != 0) return;

            const int ic = i >> 1;
            const int jc = j >> 1;

            // Bilinear weights: tx, ty in [0, 0.5]
            const Real tx = ((i & 1) == 0) ? Real(0.25) : Real(0.75);
            const Real ty = ((j & 1) == 0) ? Real(0.25) : Real(0.75);

            const Real wx0 = Real(1.0) - tx;
            const Real wx1 = tx;
            const Real wy0 = Real(1.0) - ty;
            const Real wy1 = ty;

            Real v = Real(0.0);
            v += wx0 * wy0 * a_crse(ic,   jc,   0);
            v += wx1 * wy0 * a_crse(ic+1, jc,   0);
            v += wx0 * wy1 * a_crse(ic,   jc+1, 0);
            v += wx1 * wy1 * a_crse(ic+1, jc+1, 0);

            a_fine(i,j,0) = v;
        });
    }

    child.FillBoundary(child_geom.periodicity());
}

// ============================================================================
// Spatial interpolation from parent (L0) to child (L1)
// ============================================================================
// Purpose: Prolong state variables (eta, phi) from coarse to fine level
// Usage: Initialization and buffer zone setup
//
// Key features:
// - Spatial interpolation ONLY (single time level)
// - Uses single-box coarse patch pattern (MPI-safe, multi-box safe)
// - nz=1 compatible (2D physics in 3D build)
// - Strategy: ParallelCopy gathers coarse data → interpolate → FillBoundary
//
// Related function:
// - crescent::nested::rk4::interpolateCoreDerivativeToL1Stage() handles spatial
//   + TEMPORAL interpolation for RK4 derivatives during timestepping.
//   Uses same single-box coarse patch pattern for multi-box safety.
//
// Why separate functions:
// - This function: Spatial only, used at initialization
// - CRESCENT function: Spatial + temporal blending for RK4 substages
//   Formula: D(t) = (1-c)·D^n + c·D^{n+1} where c is RK4 coefficient
// ============================================================================
void interpolateParentToChild(
    const MultiFab& parent_eta,
    const MultiFab& parent_phi,
    MultiFab&       child_eta,
    MultiFab&       child_phi,
    const Geometry& parent_geom,
    const Geometry& child_geom,
    const IntVect&  ref_ratio,
    MultiFab&       coarse_patch_scratch,
    int interpolation_order)  // 1 = bilinear, 2 = quadratic
{
    if (interpolation_order == 1) {
        prolongBilinearRef2XY(
            parent_eta,
            child_eta,
            parent_geom,
            child_geom,
            ref_ratio,
            coarse_patch_scratch);
        prolongBilinearRef2XY(
            parent_phi,
            child_phi,
            parent_geom,
            child_geom,
            ref_ratio,
            coarse_patch_scratch);
    } else if (interpolation_order == 2) {
        const int ng_fine = 2; // fill fine ghosts
        const int ng_crse = 2; // safe for quadratic (>=1 works)
        prolongQuadraticRef2XY(
            parent_eta,
            child_eta,
            parent_geom,
            child_geom,
            ref_ratio,
            ng_fine,
            ng_crse,
            coarse_patch_scratch);
        prolongQuadraticRef2XY(
            parent_phi,
            child_phi,
            parent_geom,
            child_geom,
            ref_ratio,
            ng_fine,
            ng_crse,
            coarse_patch_scratch);
    } else {
        amrex::Abort("interpolateParentToChild supports interpolation_order = 1 (bilinear) or 2 (quadratic)");
    }
}

void interpolateScalarParentToChild(
    const MultiFab& parent,
    MultiFab&       child,
    const Geometry& parent_geom,
    const Geometry& child_geom,
    const IntVect&  ref_ratio,
    MultiFab&       coarse_patch_scratch,
    int interpolation_order)
{
    if (interpolation_order == 1) {
        prolongBilinearRef2XY(
            parent,
            child,
            parent_geom,
            child_geom,
            ref_ratio,
            coarse_patch_scratch);
    } else if (interpolation_order == 2) {
        const int ng_fine = 0; // fill valid region only; caller controls child ghosts
        const int ng_crse = 2; // 3x3 coarse stencil support for quadratic interpolation
        prolongQuadraticRef2XY(
            parent,
            child,
            parent_geom,
            child_geom,
            ref_ratio,
            ng_fine,
            ng_crse,
            coarse_patch_scratch);
    } else {
        amrex::Abort("interpolateScalarParentToChild supports interpolation_order = 1 (bilinear) or 2 (quadratic)");
    }
}

IntVect getRefinementRatio(
    int lev,
    const Vector<IntVect>& ref_ratios)
{
    if (lev == 0 || lev > static_cast<int>(ref_ratios.size())) {
        return IntVect(1, 1, 1);
    }
    return ref_ratios[lev - 1];
}

} // namespace crescent::interpolation

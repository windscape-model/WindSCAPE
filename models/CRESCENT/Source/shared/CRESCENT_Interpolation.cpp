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

// Floor division for possibly negative fine indices (ghost regions).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarseIndexOf(int i, int r) noexcept
{
    return (i >= 0) ? (i / r) : (-((-i - 1) / r) - 1);
}

// Fine-cell pass of the quadratic prolongation. The Lagrange weights depend
// only on the sub-cell index m = i - r*ic of a fine cell, so they are formed
// once per call (by the same expressions as a per-cell evaluation); R > 0
// fixes the ratio at compile time, R = 0 takes it from r_runtime.
template <int R>
void prolongQuadraticFineCells(
    const MultiFab& crse_patch,
    MultiFab&       child,
    const Geometry& child_geom,
    int             r_runtime,
    int             ng_fine)
{
    constexpr int MaxRatio = 16;
    const int r = (R > 0) ? R : r_runtime;
    AMREX_ALWAYS_ASSERT(r <= MaxRatio);
    const Real inv_r = Real(1.0) / Real(r);

    // Sub-cell center s = (m + 1/2)/r - 1/2 coarse units from the coarse
    // cell center; |s| < 1/2, so the centered 3-point stencil applies.
    amrex::GpuArray<Real, 3 * MaxRatio> w{};
    for (int m = 0; m < r; ++m) {
        const Real s = (Real(m) + Real(0.5)) * inv_r - Real(0.5);
        quadraticLagrangeWeights(s, w[3 * m], w[3 * m + 1], w[3 * m + 2]);
    }

    for (MFIter mfi(child, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Region to fill: valid + ng_fine ghost (clipped to domain)
        Box bx = mfi.growntilebox(IntVect(ng_fine, ng_fine, 0));
        bx &= child_geom.Domain();

        auto const& a_fine = child.array(mfi);
        auto const& a_crse = crse_patch.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // nz=1 assumed
            if (k != 0) return;

            const int ic = coarseIndexOf(i, r);
            const int jc = coarseIndexOf(j, r);
            const int mx = i - r * ic;
            const int my = j - r * jc;
            const Real wxm1 = w[3 * mx];
            const Real wx0  = w[3 * mx + 1];
            const Real wxp1 = w[3 * mx + 2];
            const Real wym1 = w[3 * my];
            const Real wy0  = w[3 * my + 1];
            const Real wyp1 = w[3 * my + 2];

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
}

// Prolong from coarse to fine using 2D quadratic interpolation in x/y only.
// - Supports any equal x/y refinement ratio >= 2 with rz == 1 and nz = 1.
//   At r = 2 the sub-cell centers are s = -/+ 1/4.
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
    AMREX_ALWAYS_ASSERT(ref_ratio[0] >= 2 && ref_ratio[1] == ref_ratio[0]
                        && ref_ratio[2] == 1);
    AMREX_ALWAYS_ASSERT(parent.ixType() == IndexType::TheCellType());
    AMREX_ALWAYS_ASSERT(child.ixType()  == IndexType::TheCellType());
    AMREX_ALWAYS_ASSERT(child.nGrow() >= ng_fine);

    const int r = ref_ratio[0];

    // Keep the coarse valid region non-overlapping and let ParallelCopy fill
    // the coarse ghost cells directly from the parent.
    fillCoarseInterpolationPatch(parent, crse_patch, parent_geom, ng_crse);

    if (r == 2) {
        prolongQuadraticFineCells<2>(crse_patch, child, child_geom, r, ng_fine);
    } else if (r == 4) {
        prolongQuadraticFineCells<4>(crse_patch, child, child_geom, r, ng_fine);
    } else {
        prolongQuadraticFineCells<0>(crse_patch, child, child_geom, r, ng_fine);
    }

    // Fill fine ghosts across boxes/ranks (periodic)
    child.FillBoundary(child_geom.periodicity());
}

// Bilinear interpolation (1st order, simpler, works without fine ghosts).
// The cell-centered stencil supports any equal x/y integer ratio;
// the ratio-two weights remain exactly 3/4 and 1/4.
void prolongBilinearRef2XY(
    const MultiFab& parent,
    MultiFab&       child,
    const Geometry& parent_geom,
    const Geometry& child_geom,
    const IntVect&  ref_ratio,
    MultiFab&       crse_patch)
{
    AMREX_ALWAYS_ASSERT(ref_ratio[0] >= 2 &&
                        ref_ratio[1] == ref_ratio[0] &&
                        ref_ratio[2] == 1);

    const int r = ref_ratio[0];
    const Real inv_r = Real(1.0) / Real(r);

    // Bilinear only needs valid region access (no fine ghosts required)
    const int ng_crse = 1; // minimal coarse ghost for periodic BC

    fillCoarseInterpolationPatch(parent, crse_patch, parent_geom, ng_crse);

    // Bilinear interpolation kernel
    for (MFIter mfi(child, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& a_fine = child.array(mfi);
        auto const& a_crse = crse_patch.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k != 0) return;

            const int ic = coarseIndexOf(i, r);
            const int jc = coarseIndexOf(j, r);
            const Real sx =
                (Real(i - r * ic) + Real(0.5)) * inv_r - Real(0.5);
            const Real sy =
                (Real(j - r * jc) + Real(0.5)) * inv_r - Real(0.5);

            // A fine center is at s < 0 (lower children) or s > 0 (upper
            // children) relative to its containing coarse-cell center, so
            // the stencil straddles (ic-1,ic) or (ic,ic+1).
            const int ix0 = sx < Real(0.0) ? ic - 1 : ic;
            const int iy0 = sy < Real(0.0) ? jc - 1 : jc;
            const Real tx = sx < Real(0.0) ? sx + Real(1.0) : sx;
            const Real ty = sy < Real(0.0) ? sy + Real(1.0) : sy;

            const Real wx0 = Real(1.0) - tx;
            const Real wx1 = tx;
            const Real wy0 = Real(1.0) - ty;
            const Real wy1 = ty;

            Real v = Real(0.0);
            v += wx0 * wy0 * a_crse(ix0,   iy0,   0);
            v += wx1 * wy0 * a_crse(ix0+1, iy0,   0);
            v += wx0 * wy1 * a_crse(ix0,   iy0+1, 0);
            v += wx1 * wy1 * a_crse(ix0+1, iy0+1, 0);

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

} // namespace crescent::interpolation

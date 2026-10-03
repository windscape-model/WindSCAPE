#ifndef CRESCENT_NESTED_MASKS_H
#define CRESCENT_NESTED_MASKS_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include "CRESCENT_NestedDomains.H"

/**
 * @file CRESCENT_NestedMasks.H
 * @brief Blending-field generation for the active nested CRESCENT path
 *
 * The production nested algorithm uses only the derivative-blending field σ.
 * Damping-mask state is not part of the active derivative-blending workflow.
 */

namespace crescent::nested::masks {

/**
 * @brief Generate blending function σ for CRESCENT algorithm
 *
 * Purpose: Blend spectral operator results with L0 interpolation to
 *          recover accuracy at boundaries
 *
 * Definition:
 *   σ(x,y) = 0     in Ω_A (accurate domain - pure L1)
 *   σ(x,y) ∈ [0,1] in Ω_M (mixing domain - smooth blend)
 *   σ(x,y) = 1     in Ω_E (exterior - pure L0)
 *
 * Applied to derivatives/DNO outputs via blending operator:
 *   B[native_l1, core_on_l1] = (1-σ)·native_l1 + σ·core_on_l1
 *
 * Result:
 *   σ=0 in Ω_A → B = native_l1   (pure L1 spectral result)
 *   σ=1 in Ω_E → B = core_on_l1  (pure L0 interpolation)
 *
 * @param sigma_blend Output MultiFab for σ mask (single component)
 * @param geom        Geometry for L1 level
 * @param domain      CRESCENT domain structure with Ω_A, Ω_S boundaries
 * @param steepness   Reserved shape argument; this builder uses a fixed septic smoothstep profile
 */
void generateBlendingFunction(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain,
    amrex::Real steepness = 2.0
);

/**
 * @brief Generate blending function using sequential 1D approach (alternative)
 *
 * Instead of 2D Euclidean distance in corners, uses sequential 1D blending:
 * 1. First pass: blend in x-direction for all y (max of x-based σ and current σ)
 * 2. Second pass: blend in y-direction for all x (max of y-based σ and current σ)
 *
 * This preserves y-uniformity for 1D waves varying only in x.
 *
 * @param sigma_blend Output blending function MultiFab
 * @param geom        Geometry object
 * @param domain      Domain structure with Ω_A and Ω_S boundaries
 * @param steepness   Reserved shape argument; this builder uses a fixed septic smoothstep profile
 */
void generateBlendingFunctionSequential(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain,
    amrex::Real steepness = 2.0
);

/**
 * @brief Generate blending function using a 2D harmonic solve on Ω_M
 *
 * Enforces:
 *   σ = 0 on and inside Ω_A
 *   σ = 1 on the Ω_S boundary and throughout Ω_E
 *   Δσ = 0 in the interior mixing region
 *
 * This reduces ad hoc corner behavior by defining σ as a smooth 2D field.
 *
 * @param sigma_blend Output blending function MultiFab
 * @param geom        Geometry object
 * @param domain      Domain structure with Ω_A and Ω_S boundaries
 * @param steepness   Reserved shape argument
 */
void generateBlendingFunctionHarmonic(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain,
    amrex::Real steepness = 2.0
);

/**
 * @brief Generate blending function using a harmonic solve plus endpoint smoothing
 *
 * Solves the same 2D boundary-value problem as generateBlendingFunctionHarmonic(),
 * then remaps the harmonic coordinate through a septic smoothstep. This preserves
 * σ=0 on Ω_A and σ=1 on Ω_S/Ω_E while flattening derivatives near both ends.
 *
 * @param sigma_blend Output blending function MultiFab
 * @param geom        Geometry object
 * @param domain      Domain structure with Ω_A and Ω_S boundaries
 * @param steepness   Reserved shape argument
 */
void generateBlendingFunctionHarmonicSmooth(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain,
    amrex::Real steepness = 2.0
);

/**
 * @brief Generate blending function using tensor-product ramp composition
 *
 * Builds smooth 1D ramp fields in x and y, then combines them as
 *   σ = 1 - (1-σ_x)^α (1-σ_y)^α
 *
 * This always enforces:
 *   σ = 0 in Ω_A
 *   σ = 1 on and outside Ω_S
 *
 * For α=1, this gives the standard tensor-product blend. For α>1,
 * the profile stays closer to zero near Ω_A and turns on later in both
 * the side ramps and the corner regions.
 * @param sigma_blend Output blending function MultiFab
 * @param geom        Geometry object
 * @param domain      Domain structure with Ω_A and Ω_S boundaries
 * @param alpha       Corner-strength parameter (α=1 gives the standard tensor-product form)
 */
void generateBlendingFunctionTensorProduct(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain,
    amrex::Real alpha = 1.0
);

/**
 * @brief Generate blending function using a superellipse-like ramp-space radius
 *
 * Builds normalized transition coordinates t_x and t_y in [0,1] between Ω_A and Ω_S,
 * optionally remaps them as
 *   u_x = t_x^gamma,  u_y = t_y^gamma
 * combines them as
 *   r = (u_x^p + u_y^p)^(1/p)
 * and maps the clamped radius through the septic smoothstep.
 *
 * gamma=1 preserves the original 1D septic side ramps. Larger gamma delays the
 * blend turn-on near Ω_A along the side buffers. The exponent p controls how
 * rectangle-like the 2D contours are in the corner regions.
 *
 * @param sigma_blend Output blending function MultiFab
 * @param geom        Geometry object
 * @param domain      Domain structure with Ω_A and Ω_S boundaries
 * @param gamma       Side-ramp shaping exponent (1 = original side ramps)
 * @param p           Superellipse exponent (larger values are more rectangle-like)
 */
void generateBlendingFunctionSuperellipse(
    amrex::MultiFab& sigma_blend,
    const amrex::Geometry& geom,
    const domains::DomainStructure& domain,
    amrex::Real gamma = 1.0,
    amrex::Real p = 4.0
);

} // namespace crescent::nested::masks

#endif // CRESCENT_NESTED_MASKS_H

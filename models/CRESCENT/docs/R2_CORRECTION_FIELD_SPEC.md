# R2 Correction-Field Coupling (l1_coupling_mode = 2) — v1.5 Specification

**Status**: verified spec (v1.0 was adversarially reviewed by four independent
lenses on 2026-07-03; all fatal/must-fix findings are incorporated below).
**Qualification boundary (2026-08-20)**: this document specifies mode 2, but
mode 2 is no longer the production default. With the corrected JONSWAP
initializer, order-3 NS and DS runs become catastrophic at steps 607 and 507.
Three targeted candidates---additive cubic mixed terms, a compact cubic
recurrence, and a child-only order-2 DNO---remained finite to step 700 but did
not recover acceptable correction magnitude and energy fidelity. CCF is
therefore retained as an experimental mode, qualified by the monochromatic and
structured crescent cases; PTF (mode 1) is the production default.
**v1.2 (2026-07-10)**: the v1.1 background tendency `I(F_L0)` makes delta == 0
an exact invariant (the child reproduces the interpolated parent bit-for-bit);
the corrected child-assembled background tendency of §1a is now the default
(`crescent.ccf_background_tendency = 1`) and is what every published CCF
result uses. See §1a.
**v1.3 (2026-08-19)**: the corrected M=1 branch retapers the nonlocal
correction vertical velocity, honors the selected background tendency, and applies the child
pressure consistently. It is controlled during qualification by
`crescent.ccf_linear_consistency` (corrected value 1 is the default after A/B
qualification; value 0 retains the legacy branch).
**v1.4 (2026-08-20)**: the nonlinear branch applies the Dommermuth adjustment
weight consistently to the background and correction remainders
(`crescent.ccf_adjustment_consistency = 1`, default), applied pressure is kept
outside spectral projection on every path
(`crescent.consistent_applied_pressure = 1`, default), and the notation below
distinguishes the surface vertical-velocity map W from the conventional DNO G.
**v1.5 (2026-08-21)**: every background-containing metric product is
re-tapered before a child projection
(`crescent.ccf_transform_admissibility = 1`, default), and the potentials in
the W1 B terms receive the same coefficient projection as the production HOS
recursion (`crescent.ccf_w1_projection_consistency = 1`, default). The A terms
retain the unprojected potentials used to construct the first potential
coefficient.
**Status note (2026-09-05, 2026-10-03)**: the method projects the stage correction
tendency `chi*(F~0 - I(F0)) + DeltaF` and the post-step correction onto the
child modes with `max(|kx|,|ky|)` above the parent's retained cutoff
(`int(nx0/3) * dk_x0`; `applyAddedBandProjector` in
`CRESCENT_NestedCorrectionField.cpp`), so that the child's shared band evolves
exactly as the interpolated parent. Since 2026-10-03 this projection is part
of the method (`crescent.ccf_band_restrict` is a retired option). The
equations below do not show the projection; apart from it they are
unchanged.
**Status note (2026-10-02)**: the code review of 2026-10-02 retired the
compatibility switches that this document names (`ccf_background_tendency`,
`ccf_linear_consistency`, `ccf_adjustment_consistency`,
`consistent_applied_pressure`, `ccf_transform_admissibility`,
`ccf_w1_projection_consistency`, `l1_final_dno_projection`,
`l1_square_dealias`, the mask options `blending_function_mode`,
`blending_tensor_alpha` and `use_sequential_blending`, `use_c3_repair`, and the
post-timestep low-pass filter). The solver now runs only the corrected branch,
which was the default; a deck can still state that value, and any other value
stops the run. The statements below about the legacy values record the history
of this specification.
**Goal**: evolve the child as background + compact correction; no C3 repair
anywhere in the mode-2 path.

## 0. Periodicity invariant (corrected)

Every child FFT input must be either (a) **compactly supported** on the child
box, or (b) an **exact finite Fourier series** (an IFFT output, or a product
of such with a compact field). Compact support is DESTROYED by any nonlocal
spectral multiplier (|k| is Riesz-type, kernel tail ~ r^-3): vertical-velocity outputs of
compact inputs fill the box. They remain exact Fourier series (safe to
re-FFT), but their PRODUCTS with non-periodic background fields are neither
(a) nor (b) and must be re-tapered by chi before any further FFT.

Error budget: chi (septic smoothstep, exactly clamped) is C^3, so chi*delta
has Fourier coefficients that decay as ~|k|^-5 when its first unmatched trace
is the fourth derivative. A first spectral derivative then converges as the
third power of the retained band at the joins and as the fourth power at a
fixed nonzero distance from them. Compactness gives exact FFT periodicity, not
spectral convergence of arbitrary order.

## 1. Decomposition, baseline, prognostic

`I` = parent-to-child interpolation. Mode 2 REQUIRES
`state_interpolation_order == derivative_interpolation_order` (validated at
startup): one interpolation operator everywhere, or the delta telescoping
acquires a per-stage (I_state - I_deriv)(F) residue.

Child prognostic = TOTAL state u (existing RK4/update machinery).

Per-stage background baseline (delta -> 0 limit must be exact): the parent
filters its live state after every substage, so the baseline is RK-
reconstructed from PRE-filter quantities, exact by linearity of I:

```
u0_0 = I(u_saved_L0)                      (parent eta_saved/phi_saved, retained)
u0_i = I(u_saved_L0) + c_i dt I(F_L0_{i-1})   i = 1..3,  c = {0, 1/2, 1/2, 1}
u0_end = I(u_L0^{n+1})                    (parent final state, no post-filter
                                           by default; if the optional post-
                                           timestep filter is enabled, mode 2
                                           aborts at startup)
```

Correction at stage entry: `delta_i = u - u0_i`; tapered pair
`de = chi * delta_eta_i`, `dp = chi * delta_phi_i`, `chi = 1 - sigma`
(REQUIRES an exactly-clamped septic mask: `blending_function_mode == 3`
with `blending_tensor_alpha == 1`, the C^3 tensor-product mask and the
default since 2026-08-16, or the historical `blending_function_mode == 0`
Euclidean-corner mask, which is only C^1 across the four side/corner seams;
`use_sequential_blending == 0`; validated at startup).

The child two-thirds cutoffs (A1/A2 inner products, the W1 B-term potentials,
the metric products, and the post-step filter)
are per-axis square as on the parent since 2026-08-16
(`crescent.l1_square_dealias = 1`; the historical isotropic |k| cutoff is
value 0). The post-step filter's roll-off exp(-alpha (|k|/k_N)^p) is
`crescent.l1_rolloff_alpha/power` (33000/20; Hou-Li 2007 would be 36/36).

Stage tendency (written to RK component i), with the background tendency
Fbg selected by `crescent.ccf_background_tendency` (see §1a):

```
RHS_i = Fbg_i + DeltaF_i
DeltaF_i = F(u0 + chi*delta_i) - F(u0)     (tapered by construction; the
                                            mode-2 tendency deviates from the
                                            true F(u) by O((1-chi)*delta) —
                                            the tendency-level sponge, exact
                                            where chi = 1)

Fbg_i = Ftilde0_i   (ccf_background_tendency = 1, DEFAULT: §1a)
Fbg_i = I(F_L0_i)   (ccf_background_tendency = 0, legacy: delta stays 0)
```

## 1a. Background tendency (v1.2 correction; default since 2026-07-10)

With `Fbg = I(F_L0)` (v1.1), the baseline u0 advances by the SAME I(F_L0)
as the total, so the correction source `F_c(u0) - I(F_L0)` cancels
identically. DeltaF is homogeneous of degree >= 1 in delta, the child is
initialized by interpolation from L0 (delta_0 = 0), and therefore
delta == 0 for all time: the v1.1 child equals the quadratic prolongation
of the parent to round-off and evolves no sub-parent-grid content
(verified on archived runs, a fresh run, and analytically; found by the
2026-07-10 manuscript audit).

The corrected background tendency is assembled ON THE CHILD, pointwise,
from quantities already interpolated for §3 (no new captures, no child
FFT of non-periodic data — the periodicity invariant of §0 is untouched
because Fbg enters only the physical-space RK update):

```
Ftilde0_eta = - grad(psi0) . grad(eta0) + (1 + m0_d) w0
Ftilde0_psi = - g eta0 - 1/2 |grad(psi0)|^2 + 1/2 (1 + m0_d) w0sq_d
```

(pressure-free; with `consistent_applied_pressure = 1`, the default, the
applied pressure enters through -p1 after spectral projection in the assembled
RHS). The resolution commutator `C = Ftilde0 - I(F_L0)` — interpolation
does not commute with the nonlinear products — is the source that seeds
delta. On a transversely uniform parent every ingredient of Ftilde0 is
y-uniform, so the corrected background preserves transverse symmetry;
the residual child artifact arises only from the global transforms acting
on the 2D taper (measured ~1.3e-5 m on the IM case, several times below
PTF). Implementation: `runL1CorrectionFieldStage`, step 12, in
`Source/nested/CRESCENT_NestedCorrectionField.cpp`.

Post-RK4 (replaces the mode-0/1 state blend; NO C3, NO zero-mean anywhere):

```
u <- u0_end + HL[ chi * (u - u0_end) ]
```

single chi (not chi*HL[chi*.], which double-damps relative to modes 0/1);
HL = existing combined Hou-Li/dealias/low-pass filter, mean-preserving,
input compact (=> periodic). The filtered correction is a Fourier series
whose boundary leakage is O(k^-5)-small; it is the sponge + grid-scale drain.

Committed background error (default, Fbg = Ftilde0): the background
products m0_d and w0sq_d are parent-dealiased before interpolation, so
background-background content above the parent 2/3 cutoff is never
regenerated, even in Omega_A — same class as mode 1's coupling error;
document, don't hide. (Under the legacy Fbg = I(F_L0) the additional
cancellation removes the correction source entirely; see §1a.)

## 2. Parent-side per-stage captures (observational; parent untouched)

At the end of each parent RK stage: copy `{eta, phi, eta_x, eta_y, phi_x,
phi_y, phi_z}` (the post-filter state/bundle the parent's RHS actually used), and
compute ON THE PARENT GRID (spectral ops on copies):

```
w_lin   = IFFT(|k| phi_hat)          (unprojected phi, as in the phi^(1) recurrence)
lap_phi = lap(Cut23(phi))            (projected phi, as in W1 accumulation)
m0_d    = Dealias[eta_x^2 + eta_y^2] (parent-dealiased background metric)
w0sq_d  = Dealias[phi_z^2]           (parent-dealiased background w^2)
```

(The four Runge--Kutta stages and the final output state each retain 11 coarse
fields, approximately 29 MB in total at 256^2. The captured `phi` is workspace
for `w_lin` and `lap_phi`; the other 10 quantities are interpolated to the
child.)
F_L0_i needs no capture (parent RHS components are retained), nor do
eta_saved/phi_saved (retained), nor pressure (not stage-dependent).

## 3. DeltaF assembly (M = perturbation_order > 1; see §5 for M == 1)

All bg fields are stage captures interpolated with the single operator I.
Child spectral derivatives of the tapered pair are `de_x`, `de_y`, `dp_x`,
and `dp_y`. Define
`dm_raw = 2(eta0_x de_x + eta0_y de_y) + de_x^2 + de_y^2`,
`dm = Dealias_child[chi * dm_raw]`, and `dwt = chi * dw` (§4):

```
DeltaF_eta = - (psi0_x de_x + psi0_y de_y) - (dp_x eta0_x + dp_y eta0_y)
             - (dp_x de_x + dp_y de_y)
             + (1 + I(m0_d) + dm) * dwt + dm * w0

DeltaF_psi = - g de - (psi0_x dp_x + psi0_y dp_y) - 0.5 (dp_x^2 + dp_y^2)
             + 0.5 (1 + I(m0_d)) * Dealias_child[2 w0 dwt + dwt^2]
             + 0.5 * dm_w2
  where dm_w2 = Dealias_child[ chi * dm_raw
                               * (I(w0sq_d) + 2 w0 dwt + dwt^2) ]
```

The explicit `chi` factors preceding these projections restore compact support
after a transformed correction derivative has been multiplied by a
non-periodic background field. The projected outputs may fill the box, but
they are exact child Fourier series until the final taper below.

Final per-stage dealias (parity with the parent's rhs_eta/rhs_phi points):

```
DeltaF_eta <- Dealias_child[chi * DeltaF_eta]
DeltaF_psi <- Dealias_child[chi * DeltaF_psi]   then
DeltaF_psi <- DeltaF_psi - (p1 - I(p0))         (pressure delta: added AFTER
                                                 the dealias, never through a
                                                 child FFT; zero in all
                                                 current cases)
```

Dealiasing here is parity with the parent recipe (metric, w^2, rhs_eta,
rhs_psi points), not exactness: compactness licenses the FFTs; the 2/3 rule
gives the same aliasing-control standard as the parent. (Optional
improvement unavailable to the parent: pre-truncate de/dp to the child 2/3
band before products.)

## 4. Correction vertical velocity dw

```
dw = W_M[de](dp)     (existing child vertical-velocity machinery on the tapered pair;
                      internals audited: eta-powers of de are compact,
                      phi-expansion FFT inputs are compact-or-Fourier-series,
                      per-order 2/3 cutoffs as in the parent chain)
   + A1 + A2 + B1 + B2      (order-2 bg-delta cross terms, only if M >= 2):
     A1 = - IFFT( |k| Cut23[ FFT( chi * ( eta0 * IFFT(|k| dp_hat) ) ) ] )
     A2 = - IFFT( |k| Cut23[ FFT( de * I(w_lin0) ) ] )
     B1 = - eta0 * lap(Cut23(dp))
     B2 = - de * I(lap(Cut23_0(phi0)))
then dwt = chi * Cut23[chi * dw]    (taper -> final complete-sum projection
                                     -> taper; safe before every background
                                     product)
```

Here `W[eta] psi = phi_z|z=eta` is the surface vertical-velocity map. The
conventional DNO is reconstructed as
`G[eta] psi = (1 + |grad eta|^2) W[eta] psi - grad eta . grad psi`.

Notes: A1's inner product is chi-tapered (eta0 is non-periodic and
IFFT(|k| dp_hat) is box-filling — the product needs re-compactification);
A1/A2 inner FFTs receive the 2/3 cutoff BEFORE the outer |k| (mirroring the
parent's per-order cutoff, which cannot remove content already aliased into
the retained band). The A terms use the unprojected surface potentials because
they construct the first potential coefficient; the B terms use the projected
zeroth coefficients that enter the vertical-velocity accumulation.
Disjointness: W_M[de](dp) is all-delta; A/B are strictly
one-bg-one-delta; bg-bg lives in w0 = I(phi_z capture). v1 omits cubic mixed
terms (bg-bg-delta, bg-delta-delta): O((k a0)^2) and O(k a0 kd) relative on
dw (~1% at ak ~ 0.1).

## 4a. Dommermuth nonlinear adjustment

For `adjustment_periods > 0`, let `Ta` be that many peak periods and
`A(t) = 1 - exp[-(t/Ta)^4]` while `t < 2 Ta` (the implementation returns
exactly one thereafter). With `ccf_adjustment_consistency = 1`, the default,
the child uses the same linear-plus-ramped-remainder decomposition as L0:

```
Fbg_adjusted   = L(u0)     + A(t) [Ftilde0 - L(u0)]
DeltaF_adjusted = L(delta) + A(t) [DeltaF  - L(delta)]
RHS             = Fbg_adjusted + DeltaF_adjusted - p1
```

The compact linear correction follows the same taper, final vertical-velocity
projection (`l1_final_dno_projection`), and dealias chain as the nonlinear
increment. Applied pressure is subtracted
after the ramp and is never scaled. Setting `ccf_adjustment_consistency = 0`
with a positive adjustment time is rejected at startup rather than silently
mixing ramped parent physics with an unramped correction.

## 5. M == 1 (linear) branch

The parent's M=1 branch is strictly linear. The consistent mode-2 branch uses:

```
DeltaF_eta = chi * IFFT(|k| dp_hat)
DeltaF_psi = - g de - (p1 - I(p0))
```

The corrected background selection is honored explicitly: the default uses
`Ftilde0 = (I(w0), -g I(eta0) - p1)`, while the legacy background uses
`I(F_L0)` plus the displayed pressure increment. There are no cross terms, no
nonlinear table, and no dealias calls. During A/B qualification, setting
`ccf_linear_consistency = 0` preserves the pre-v1.3 branch.

## 6. Mode-2 code-path requirements (from the review)

1. The per-stage periodic-scratch preparation (C3 + Hou-Li on the total
   state) is NOT invoked in mode 2 — the stage branch must skip it entirely,
   not merely ignore its outputs (with use_c3_repair=0 it would FFT a
   non-periodic total state).
2. The post-RK4 path in mode 2 is the §1 taper — it must NOT call the
   mode-0/1 preparePeriodicL1StateScratch (which applies C3 unconditionally,
   ignoring use_c3_repair).
3. Outputs/energy: the final-output blended-derivative recompute (C3-based)
   is replaced in mode 2 by the split: output bundle =
   I(parent np1 bundle) + spectral derivatives of the FINAL tapered
   correction (+ dw_end), written into the same l1_blended_derivatives
   storage so tecplot/energy/edge diagnostics work unchanged. The standard-
   bundle recompute (periodic surface-map evaluation on the raw total state) is skipped for the
   child in mode 2 for the same reason. Since 2026-08-16 the mode-2 output
   dw_end includes the order-2 cross terms A1+A2+B1+B2 formed with the
   FINAL parent background (u_L0^{n+1} and its bundle) and the same chi
   re-taper as the stages, so the output phi_z is the stage-consistent one
   (before, it held only W_M[de]dp: an energy/diagnostic-only defect). Energy
   reconstructs the conventional DNO pointwise from this bundle as
   `G psi = (1 + |grad eta|^2) w - grad eta . grad psi`.
4. Seam diagnostics in mode 2 monitor the fields that actually enter FFTs
   (de, dp, chi*DeltaF, the A1 inner product), not the unused periodic
   scratch.
5. Startup validation (abort otherwise): equal interpolation orders;
   septic-clamped mask (blending_function_mode==3 with tensor alpha 1, or 0;
   use_sequential_blending==0);
   post-timestep low-pass filter disabled; use_c3_repair may be 0 or 1 (it
   is simply never consulted in the mode-2 path).

## 7. Verification obligations

1. delta -> 0: with the §1 pre-filter RK-reconstructed baseline, the child
   tracks I(parent) exactly through all stages and the post-step (the
   substage-filter mismatch of v1.0 is eliminated by construction).
2. chi == 1, equal resolutions: total dynamics = full HOS to the v1
   vertical-velocity-map cutoff and the captured-dealias parity of §3.
3. With M == 2 and the taper suppressed, the projected W1 construction equals
   the complete discrete vertical-velocity increment to round-off.
4. M == 1 reduces to the parent's linear system exactly.
5. Every FFT input satisfies §0; enforced by construction and monitored for
   the principal transformed fields by the §6.4 diagnostics.
6. Omega_A mean watch: mean(DeltaF_eta) is generically nonzero and nothing
   pins the Omega_A mean of delta_eta; monitor for secular drift in the
   validation runs (add mean(delta_eta over Omega_A) to the oscillation CSV
   in mode 2).

## 8. Acceptance

Mode 2 is checked against mode 1 on the retained monochromatic and structured
crescent cases (its structural advantages remain no C3 repair and no
total-state splice); B1b seam/artifact metrics should improve or hold. It does
not satisfy the corrected order-3 NS/DS broadband acceptance cases and must not
be presented as qualified there. The production regression gate (default =
mode 1) must remain bit-for-bit.

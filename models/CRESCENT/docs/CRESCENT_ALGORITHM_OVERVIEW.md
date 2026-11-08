# CRESCENT Algorithm Overview

*CRESCENT: A block-structured, adaptive-mesh framework for ocean surface wave modeling*

CRESCENT is a C++17/AMReX research code for high-order spectral simulation of
nonlinear ocean free-surface waves. The primary prognostic fields are the
surface elevation `eta` and the surface velocity potential `phi`. The core path
advances a single periodic domain; the nested path adds one-way L0-to-L1
coupling for localized high-resolution wave dynamics, and offers two
interchangeable coupling methods, the periodized total-field method (PTF,
the production default) and the compact correction-field method (CCF).

## Solver Structure

The active source tree separates major responsibilities into focused areas:

- `Source/core/`: single-domain free-surface evolution, spectral derivatives,
  surface vertical-velocity evaluation, and physical RHS construction.
- `Source/nested/`: nested L1 setup, the two coupling methods (correction-field
  and periodized total-field), C3-target periodization, tendency blending, nested RK4
  evolution, and nested diagnostics.
- `Source/solver/`: AMReX driver ownership, runtime input parsing, setup,
  checkpoint/restart, output, run loop, and reporting.
- `Source/shared/`: common state views, mesh views, runtime state, spectral
  workspaces, interpolation, and utility code.
- `Source/external/`: non-owning surface exchange views intended for
  coupling with external solvers.

This layout keeps the top-level `CRESCENT` class as the AMReX-facing driver while
moving numerical operations into domain-specific modules.

## Core Free-Surface Evolution

The core solver path advances the periodic free-surface state using high-order
spectral operations. A typical step computes horizontal derivatives and the
surface vertical velocity through the HOS vertical-velocity series, reconstructs
the conventional Dirichlet-to-Neumann value
`G psi = (1 + |grad eta|^2) w - grad eta . grad psi` when needed, builds the
Zakharov free-surface RHS, and advances `eta` and `phi` with the active time
integration path.

Spectral workspaces and scratch fields are allocated during setup and reused in
the timestep loop. This avoids hot-path heap allocation and keeps AMReX
`MultiFab` ownership explicit.

## Nested Refinement

The nested L1 algorithm refines a non-periodic rectangular region embedded
inside a periodic L0 wave solution. The coupling is one-way: L0 evolves
independently, and L1 reads L0 state and tendencies without feeding back to L0.
Both levels share the time step, and L0 advances first within each step so that
its stage tendencies are available to L1 exactly. Nested runs use a fixed AMR
hierarchy after initialization; the spectral and nested scratch workspaces are
allocated for the initialized level layout and are rebuilt during restart, not
during in-run regridding. Parent-to-child transfer uses a single quadratic
prolongation operator by default, and the parent-to-child refinement ratio may
be any equal integer ratio of two or more in the horizontal directions. The
optional ratio-two, order-1 transfer uses the correct cell-centered bilinear
stencil; the historical half-cell-shifted stencil was retired on 2026-10-02.

Three nested zones organize the coupling through a smooth mask `sigma`: the
accurate core `Omega_A` where `sigma = 0`, a mixing annulus `Omega_M` where
`sigma` ramps from zero to one through a C3 septic smoothstep, and an exterior
buffer `Omega_E` where `sigma = 1`. The complementary taper is `chi = 1 - sigma`.

### The periodicity discipline

Both methods are organized around one invariant: every field passed to a child
FFT must be either compactly supported on the child box or an exact finite
Fourier series on it. Non-periodic fields, which include the live total state,
every field interpolated from L0, and any `sigma`-blended combination, are
touched only pointwise in physical space. The invariant is not preserved by the
nonlocal `|k|` multipliers of the surface vertical-velocity map, whose outputs
fill the box even for compact inputs and must therefore be re-tapered by `chi`
before they multiply non-periodic background fields. Two further requirements
follow: child-side filters preserve the `k = 0` mode, because a window of a
larger sea has a nonzero mean, and every term of an assembled tendency derives
from one consistent field pair, because discrete energy stability rests on the
Hamiltonian structure of the equations.

### Compact correction-field method (CCF, mode 2)

The correction-field method evolves the deviation of the child solution from
the interpolated parent field rather than the total state. All non-periodicity
resides in the parent background, which is never transformed on the child, while
the correction is tapered by `chi` to compact support and is therefore periodic
by construction, so no periodization step is required. Per RK stage, the child
reconstructs the stage-exact interpolated parent baseline from pre-filter parent
quantities, forms the tapered correction, evaluates the correction tendency
increments (including the correction vertical-velocity map built from
the parent map, a full-order recursion on the correction pair, and explicit
cross terms), re-tapers every background-containing product before a child
projection, and assembles the stage tendency as the interpolated parent
tendency plus the projected correction tendency. The linear-in-elevation cross terms follow the parent
recursion's coefficient ordering: the outer-map terms use the unprojected
surface potentials, while the Laplacian terms use the projected zeroth
coefficients. The background tendency is assembled on the child from
the interpolated parent fields and parent-dealiased products (pointwise, never
transformed). This choice matters: with the interpolated parent tendency as
background (the behavior before 2026-07-10, retired in 2026-10), the
correction source cancels against the baseline, and the child reproduces the
interpolated parent exactly, evolving no correction at all. With the
child-assembled background, the difference between the two background
tendencies (the resolution commutator) seeds the
correction. After the RK combine, the child state is projected back onto the
decomposition, which damps and filters the correction alone. The method has no
periodization tuning parameters. Its discrete specification is given in
`docs/R2_CORRECTION_FIELD_SPEC.md`.

The correction holds only the content the parent cannot represent: the stage
correction tendency and the post-step correction are projected onto the child
modes above the parent's retained cutoff, so that the band shared with the
parent evolves exactly as the interpolated parent. A correction in the shared
band would alter the parent's resolved waves, which a one-way coupling cannot
keep consistent with the parent.

### Periodized total-field method (PTF, mode 1, production default)

The periodized total-field method evolves the total child state and restores
FFT admissibility at every RK stage. Scratch copies of `eta` and `phi` are made
C3-target periodic by a regularized least-squares repair confined to the exterior buffer
(the physical state is never modified), the copies are filtered by the
mean-preserving combined spectral filter, and the child tendency is computed
entirely from the periodized copies and dealiased exactly as on L0. The child
tendency is then blended with the stage-exact interpolated L0 tendency through
the mask, `RHS = (1 - sigma) F_L1 + sigma I(F_L0)`, and the RK4 update is applied
in physical space. After the RK stages, a post-step relaxation replaces the
child state by a `sigma`-blend of the repaired-and-filtered copy and the
interpolated L0 state, which both enforces consistency near the boundary and
provides the grid-scale dissipation of the non-periodic live state.

## C3-target Periodization (PTF only)

FFT-based differentiation requires periodic input. In the periodized total-field
method the L1 region is not physically periodic, so CRESCENT applies a
regularized repair that targets continuity of the value and first three
derivatives to scratch copies of the nested state before spectral operations. The physical L1
state is not modified by this repair step; the repaired copies are used only to
compute periodic spectral derivatives and surface vertical velocity safely, within the
exterior buffer `Omega_E`. The correction-field method requires no periodization
and does not use this path. At the baseline regularization the wrap conditions
are not exact; seam diagnostics quantify the residual, while the unregularized
exact-interface limit realizes the formal C3 decay law.

## Runtime Configuration

Input decks use the `crescent.` prefix. Important runtime groups are:

- Time control: `max_step`, `stop_time`, `crescent.dt`, `crescent.dt_override`.
- Wave initialization: `crescent.wave_type`, `crescent.wind_speed_10m`,
  `crescent.fetch`, `crescent.wave_steepness`, `crescent.wave_mode_x`.
  For JONSWAP, `wave_steepness` is inert: the spectrum fixes the amplitude,
  with `crescent.jonswap_amplitude_scale` as an explicit multiplier. The
  construction is normalized, uses signed angles, and represents both
  propagation senses of each real-FFT conjugate pair.
  `crescent.use_directional_spectrum=0` makes the spectrum long-crested in
  `+x`. Stokes angles keep all four quadrants.
- Spectral controls: `crescent.perturbation_order`, `crescent.beta1_filter`,
  `crescent.beta2_filter`, and the optional per-level override
  `crescent.beta1_filter_l1`. A positive
  `crescent.peak_wavenumber_override` replaces the initializer's filter
  peak. The `beta1_filter`/`beta2_filter` pair also
  serves as the spectral wave-breaking dissipation model of Xiao et al. (2013)
  when set to their breaking-calibrated values.
- Nested refinement: `crescent.use_crescent_algorithm`,
  `crescent.l1_coupling_mode` (`1` = PTF, default; `2` = CCF),
  `crescent.L1_box_lo`/`crescent.L1_box_hi`, `crescent.roi_accurate_fraction`,
  `crescent.roi_spectral_fraction`, and, for PTF, the
  `crescent.c3_boundary_penalty` and `crescent.c3_coefficient_regularization`
  periodization weights. The parent-to-child refinement ratio is set through
  `amr.ref_ratio_vect`. The compatibility switches of the 2026-07 and
  2026-08 fixes and `crescent.ccf_band_restrict` were retired in 2026-10: a
  deck may state their supported value, and any other value stops the run.
- Steep-sea initialization: `crescent.adjustment_periods` applies a
  Dommermuth-style time ramp of the nonlinear tendency (nonzero enables it;
  supported by the core, PTF, and corrected CCF paths).
- Output controls: `crescent.plot_file`, `crescent.plot_int`,
  `crescent.tecplot_int`, `crescent.write_omega_a_edge_diagnostics`.
- Checkpoint/restart: `crescent.check_file`, `crescent.check_int`,
  `crescent.restart`.

Representative inputs are provided in `test_case/`.

## GPU Execution

Both nested methods run on GPUs. Configure with `-DCRESCENT_ENABLE_CUDA=ON`
(plus `-DCMAKE_CUDA_ARCHITECTURES=<arch>`); the analogous
`CRESCENT_ENABLE_HIP`/`CRESCENT_ENABLE_SYCL` switches select the other AMReX
backends. Every per-cell operation of the core solver and the correction-field
path is an AMReX device kernel, transforms go through `amrex::FFT::R2C` (cuFFT
on CUDA), and the C3 interface repair of the total-field method runs as a
batched device kernel that solves one grid line's least-squares system per
thread. Output and diagnostics stage through pinned host memory. GPU and CPU
runs of the same deck agree to machine precision, not bit for bit: FFT
algorithms, parallel reductions, and fused multiply-add contraction change the
rounding.

## Checkpoint And Restart

The checkpoint path writes enough information to restart the time integration
without relying on initialization-only state. Checkpoints contain the
free-surface state fields, step/time counters, runtime reporting state, level
geometry and layout metadata, the complete wave-initialization, spectral, and
nested evolution configuration, and the metadata needed to rebuild derived
fields and scratch workspaces. The checkpoint format (version 10) restores
the checkpointed evolution settings rather than silently inheriting
conflicting values from the restart deck; the slots of retired options hold
their supported value.

Derived fields such as `phi_z`, horizontal derivatives, spectral workspaces,
and C3/nested runtime structures are rebuilt or recomputed during restart rather
than serialized as authoritative state. This keeps checkpoint files focused on
the physical solution and reproducible runtime metadata.

The restart path supports both single-domain runs and nested runs. The input
file used for restart must supply a compatible hierarchy; max step, stop time,
and output controls may be changed for the continuation.

## External Exchange Interface

The external exchange API is intentionally lightweight. `CRESCENT` implements
`crescent::external_exchange::SurfaceExchangeProvider`, which exposes non-owning
AMReX views of the current level state:

- `eta`: free-surface elevation.
- `phi`: free-surface potential.
- `pressure`: pressure field.
- `Geometry`, `BoxArray`, and `DistributionMapping`: AMReX metadata needed to
  interpret the fields.

The interface is designed for coupling infrastructure. It does not define
a coupling algorithm, schedule, interpolation policy, or feedback mechanism.

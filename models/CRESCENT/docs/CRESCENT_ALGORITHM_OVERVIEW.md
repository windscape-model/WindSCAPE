# CRESCENT Algorithm Overview

*CRESCENT: A block-structured, adaptive-mesh framework for ocean surface wave modeling*

CRESCENT is a C++17/AMReX research code for high-order spectral simulation of
nonlinear ocean free-surface waves. The primary prognostic fields are the
surface elevation `eta` and the surface velocity potential `phi`. The core path
advances a single periodic domain; the nested path adds one-way L0-to-L1
coupling for localized high-resolution wave dynamics.

## Solver Structure

The active source tree separates major responsibilities into focused areas:

- `Source/core/`: single-domain free-surface evolution, spectral derivatives,
  DNO evaluation, and physical RHS construction.
- `Source/nested/`: nested L1 setup, C3 periodization, derivative blending,
  nested RK4 evolution, and nested diagnostics.
- `Source/solver/`: AMReX driver ownership, runtime input parsing, setup,
  checkpoint/restart, output, run loop, and reporting.
- `Source/shared/`: common state views, mesh views, runtime state, spectral
  workspaces, interpolation, and utility code.
- `Source/external/`: non-owning surface exchange views intended for future
  coupling with external solvers.

This layout keeps the top-level `CRESCENT` class as the AMReX-facing driver while
moving numerical operations into domain-specific modules.

## Core Free-Surface Evolution

The core solver path advances the periodic free-surface state using high-order
spectral operations. A typical step computes horizontal derivatives and the
vertical surface velocity through the HOS Dirichlet-to-Neumann operator, builds
the Zakharov free-surface RHS, and advances `eta` and `phi` with the active time
integration path.

Spectral workspaces and scratch fields are allocated during setup and reused in
the timestep loop. This avoids hot-path heap allocation and keeps AMReX
`MultiFab` ownership explicit.

## Nested Repair-Compute-Blend Path

The nested L1 algorithm is designed for a non-periodic refined region embedded
inside a periodic L0 wave solution. The active approach is one-way coupled: L0
evolves independently, and L1 reads L0 state and derivatives without feeding
back to L0.

Nested runs use a fixed AMR hierarchy after initialization. The spectral, C3,
and nested scratch workspaces are allocated for the initialized level layout and
are rebuilt during restart, not during in-run regridding.

For each nested RK stage:

1. Save or advance the L1 RK substage state before derivative evaluation.
2. Copy L1 `eta` and `phi` into periodic scratch fields.
3. Apply C3 repair to the scratch fields so FFT operations see periodic data.
4. Apply spectral filtering only to periodic fields where filtering is valid.
5. Compute L1 spectral derivatives and vertical velocity on the repaired fields.
6. Interpolate the appropriate L0 derivative snapshot to the L1 grid.
7. Blend L1 and interpolated L0 derivatives through the sigma field.
8. Compute the nested free-surface RHS from the blended derivatives.

After the RK stages, the nested state is filtered and blended with the L0 state
in the exterior coupling zone. This keeps the accurate core of L1 governed by
the refined solver while enforcing consistency with L0 near the outer boundary.

## C3 Periodization

FFT-based differentiation requires periodic input. The L1 refined region is not
physically periodic, so CRESCENT applies C3 repair to scratch copies of the nested
state before spectral operations. The physical L1 state is not modified by this
repair step. The repaired copies are used only to compute periodic spectral
derivatives and DNO quantities safely.

The coupling zones are commonly described as:

- `Omega_A`: accurate L1 core where sigma is zero or near zero.
- `Omega_M`: transition region where L1 and interpolated L0 information blend.
- `Omega_E`: exterior region dominated by L0 information and used by the C3
  repair buffer.

## Runtime Configuration

Input decks use the `crescent.` prefix. Important runtime groups are:

- Time control: `max_step`, `stop_time`, `crescent.dt`, `crescent.dt_override`.
- Wave initialization: `crescent.wave_type`, `crescent.wind_speed_10m`,
  `crescent.fetch`, `crescent.wave_steepness`, `crescent.wave_mode_x`.
- Spectral controls: `crescent.perturbation_order`, `crescent.beta1_filter`,
  `crescent.beta2_filter`.
- Output controls: `crescent.plot_file`, `crescent.plot_int`,
  `crescent.tecplot_int`.
- Checkpoint/restart: `crescent.check_file`, `crescent.check_int`,
  `crescent.restart`.

Representative inputs are provided in `test_case/`.

## Checkpoint And Restart

The checkpoint path writes enough information to restart the time integration
without relying on initialization-only state. Checkpoints contain the
free-surface state fields, step/time counters, runtime reporting state, level
geometry and layout metadata, nested configuration, and the metadata needed to
rebuild derived fields and scratch workspaces.

Derived fields such as `phi_z`, horizontal derivatives, spectral workspaces,
and C3/nested runtime structures are rebuilt or recomputed during restart rather
than serialized as authoritative state. This keeps checkpoint files focused on
the physical solution and reproducible runtime metadata.

The restart path supports both single-domain runs and nested runs. The input
file used for restart must be compatible with the checkpointed hierarchy and
configuration.

## External Exchange Interface

The external exchange API is intentionally lightweight. `CRESCENT` implements
`crescent::external_exchange::SurfaceExchangeProvider`, which exposes non-owning
AMReX views of the current level state:

- `eta`: free-surface elevation.
- `phi`: free-surface potential.
- `pressure`: pressure field.
- `Geometry`, `BoxArray`, and `DistributionMapping`: AMReX metadata needed to
  interpret the fields.

The interface is designed for future coupling infrastructure. It does not define
a coupling algorithm, schedule, interpolation policy, or feedback mechanism.

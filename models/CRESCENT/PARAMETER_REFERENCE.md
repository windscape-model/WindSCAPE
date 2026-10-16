# CRESCENT Input Parameter Reference

Complete reference for the input parameters read by the CRESCENT wave solver
(High-Order Spectral, HOS). Generated from the parameter readers in
`Source/solver/CRESCENT_SetupParameters.cpp` and
`Source/nested/CRESCENT_NestedAlgorithmInputs.cpp`; defaults are taken from the
state structs (`CRESCENT_RuntimeControl.H`, `CRESCENT_WaveState.H`,
`CRESCENT_NestedConfig.H`). Last refreshed 2026-10-03 against the standalone
sync at `f15a84e`; updated for the CRESCENT 1.1 release tree `b34df45`
(`crescent.ccf_band_restrict` retired). Changes since the 2026-09-04 refresh (`f228bf6`):

- Retired inputs (2026-10-02 code review): 50 legacy or exploratory switches
  and the historical derivative blend (`l1_coupling_mode = 0`) were removed
  from the solver. A deck can still give a retired switch its former default;
  any other value stops the run with a message that names commit `f39cacf`.
  `l1_coupling_mode = 0` stops the run. [Retired inputs](#retired-inputs)
  lists each switch with the one value that the solver still accepts.
- `crescent.ccf_band_restrict` (added 2026-09-05) was retired on 2026-10-03: the projection onto the
  added band is part of CCF; a deck may still set 1.
- No longer read (2026-09-25): `crescent.enable_initial_c3_effect_diagnostic`
  and `crescent.enable_nested_initialization_validation`; their only effect was
  diagnostic output that was removed.

Earlier changes that still apply (2026-07-17 to 2026-09-04): every parameter is
range-checked at startup; the derived quantities `alpha`, `peak_frequency`,
`peak_wavelength` and `swell_wavenumber` are rejected as inputs; `crescent.dt`
is required for the `stokes` / `crescent` wave types and rejected for
`jonswap`, which uses `crescent.dt_override` (or the automatic T_peak/20); the
defaults of `beta1_filter` (8 -> 9), `c3_boundary_penalty` (0.001 -> 0.01),
`c3_coefficient_regularization` (1 -> 10) and `l1_coupling_mode` (2 -> 1, PTF)
changed.

---

## Parameter prefix (read this first)

CRESCENT reads its solver parameters through an AMReX `ParmParse` group named by
the **`crescent`** prefix. How the prefix applies depends on how CRESCENT is run:

| Run mode | Solver parameters | `amr.*` / `geometry.*` |
|----------|-------------------|------------------------|
| **Standalone** (`crescent` executable) | `crescent.<name>` | **unprefixed** (`amr.*`, `geometry.*`) |
| **Coupled** (`windscape`, ERF + CRESCENT) | `crescent.<name>` | **prefixed** (`crescent.amr.*`, `crescent.geometry.*`) |

A handful of top-level run-control parameters (`max_step`, `stop_time`,
`start_time`) are read from the **global** (unprefixed) `ParmParse` in both modes.
In coupled mode the run length comes from ERF; CRESCENT's own `max_step` /
`stop_time` are not consulted.

> Earlier revisions of this solver used a `hos.` prefix and a `use_stokes_waves`
> flag. **Both are obsolete.** The prefix is now `crescent.` and the wave model
> is selected with `crescent.wave_type` (a string). Parameters set under `hos.`
> are silently ignored by the current code.

---

## Parameter Categories

1. [Time Control](#time-control)
2. [Domain (AMR & Geometry)](#domain-amr--geometry)
3. [Spectral / HOS Solver](#spectral--hos-solver)
4. [Wave Initialization](#wave-initialization)
5. [Output Control](#output-control)
6. [Checkpoint / Restart](#checkpoint--restart)
7. [Nested AMR (CRESCENT algorithm)](#nested-amr-crescent-algorithm)
8. [Retired inputs](#retired-inputs)

---

## Time Control

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `max_step` | int | Maximum number of timesteps (global, unprefixed) | `INT_MAX` | `>= 0` |
| `stop_time` | real | Simulation end time, s (global, unprefixed) | `REAL_MAX` | finite, `>= start_time` |
| `start_time` | real | Simulation start time, s (global, unprefixed) | `0.0` | finite |
| `crescent.dt` | real | Fixed timestep, s, for `wave_type = stokes` or `crescent` | none (required) | positive finite; **rejected** for `jonswap` |
| `crescent.dt_override` | real | JONSWAP only: explicit timestep, s | T_peak/20 (automatic) | positive finite; **rejected** for other wave types |

The run ends at whichever of `max_step` / `stop_time` is reached first.

---

## Domain (AMR & Geometry)

Standard AMReX `AmrCore` parameters. In **standalone** runs these are unprefixed;
in **coupled** runs prepend `crescent.` (e.g. `crescent.amr.n_cell`).

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `amr.n_cell` | int int int | Grid cells (x, y, z) | Required |
| `amr.max_level` | int | Maximum AMR level: `0` = single level, `1` = one nested level (higher values are rejected) | `0` |
| `amr.ref_ratio_vect` | int int int | Refinement ratio for the nested level (use `2 2 1` or `4 4 1`) | AMReX default |
| `amr.max_grid_size` | int | Maximum grid (box) size | AMReX default |
| `amr.blocking_factor` | int | Grid blocking factor | AMReX default |
| `geometry.prob_lo` | real real real | Lower domain corner, m | `0 0 0` |
| `geometry.prob_hi` | real real real | Upper domain corner, m | none |
| `geometry.prob_extent` | real real real | Domain size, m (alternative to `prob_hi`) | none |
| `geometry.is_periodic` | int int int | Periodic BC flags (x, y, z) | Required (use `1 1 0`) |

**Notes**
- The HOS spectral method is single-level on L0; set `amr.max_level = 0` unless
  you are using the nested CRESCENT algorithm (see below).
- Use periodic in x and y for ocean waves: `geometry.is_periodic = 1 1 0`.
- `n_cell` in x/y should be powers of 2 for FFT efficiency.
- Domain size sets the wavenumber resolution dk = 2 pi / L (computed automatically).
- Nested runs (`amr.max_level = 1`) currently require exactly one MPI rank.

---

## Spectral / HOS Solver

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.perturbation_order` | int | HOS expansion order M | `3` | `>= 1` |
| `crescent.beta1_filter` | real | Spectral low-pass filter width (manuscript baseline 9 since 2026-08-16; was 8) | `9.0` | positive finite |
| `crescent.beta1_filter_l1` | real | L1 override of `beta1_filter`, so each level can damp its own cutoff-pileup band; `-1.0` = use `beta1_filter` | `-1.0` | `-1` or positive finite |
| `crescent.beta2_filter` | real | Spectral filter parameter 2 | `12.0` | positive finite |
| `crescent.l1_rolloff_alpha` | real | L1 spectral roll-off strength in exp(-alpha (\|k\|/k_N)^p) (historical 33000; Hou-Li 2007 uses 36) | `33000.0` | positive finite |
| `crescent.l1_rolloff_power` | int | L1 spectral roll-off exponent p (historical 20; Hou-Li 2007 uses 36) | `20` | `>= 1` |
| `crescent.adjustment_periods` | real | Dommermuth (2000) nonlinear ramp time T_a, in peak periods; ramps the nonlinear terms in from a linear initialization (needed for steep seas started linearly); `0` = off | `0.0` | non-negative finite |
| `crescent.peak_wavenumber_override` | real | Peak wavenumber for the filter; `-1.0` = compute automatically | `-1.0` | `-1` or positive finite |

### Perturbation Order

`perturbation_order` (M) trades accuracy for cost. Higher M captures steeper,
more nonlinear waves:

| M | Regime | Typical use |
|---|--------|-------------|
| 1 | Linear | Testing only |
| 2 | Weakly nonlinear | Mild seas |
| 3 | Nonlinear | **Recommended default** |
| 4-5 | Highly nonlinear | Steep / near-breaking waves |
| 6+ | Extreme | Research |

> Wavenumber spacing (dk_x, dk_y) is **computed internally** from the domain size;
> it is not an input parameter.

---

## Wave Initialization

The initialization model is chosen with `crescent.wave_type`.

### Common

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.wave_type` | string | `"jonswap"`, `"stokes"`, or `"crescent"` (the legacy names `"monochromatic"` / `"focused_group"` are not dispatched by the current initializer) | `"jonswap"` | one of the three |
| `crescent.wave_steepness` | real | Wave steepness ka (Stokes / crescent amplitude control; not used by JONSWAP) | `0.1` | non-negative finite |
| `crescent.initial_angle` | real | Initial/mean wave direction, rad | `0.0` | finite |
| `crescent.wave_angle` | real | Stokes propagation angle, rad (0 = +x, pi/2 = +y) | `0.0` | finite |

### Stokes waves (`crescent.wave_type = stokes`)

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.wave_mode_x` | int | Wave mode index: the carrier wavenumber is k = `wave_mode_x` * 2 pi / L_x (rotated by `wave_angle`) | `5` | `>= 1` |

Wavelength relation: lambda = 2 pi / k. Amplitude from steepness: a = ka / k.
`crescent.swell_wavenumber` is **not** an input: no initializer reads it, and the
current solver aborts when it is present.

### Crescent (horseshoe) waves (`crescent.wave_type = crescent`)

A 5th-order Stokes carrier along +x plus an oblique free-wave perturbation
(the class-II / McLean instability seed that develops into crescent waves).
The carrier wavenumber and steepness come from `wave_mode_x` / `wave_steepness`;
the ratios below are the perturbation wavenumbers in carrier units.

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.crescent_eps` | real | Perturbation amplitude relative to the carrier amplitude a | `0.16` | finite |
| `crescent.crescent_kx_ratio` | real | Perturbation k_x / carrier k (streamwise 3/2 subharmonic) | `1.5` | `> 0` |
| `crescent.crescent_ky_ratio` | real | Perturbation k_y / carrier k (transverse oblique mode) | `1.23` | `> 0` |
| `crescent.crescent_stokes_order` | int | Stokes order of the analytic carrier | `5` | `1-5` |

### JONSWAP spectrum (`crescent.wave_type = jonswap`)

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.wind_speed_10m` | real | Wind speed at 10 m, m/s | `0.0` | positive finite, required |
| `crescent.fetch` | real | Fetch length, m | `0.0` | positive finite, required |
| `crescent.jonswap_random_seed` | long | Phase seed for reproducible init; `-1` = AMReX default seeding | `12345` | `>= -1` |
| `crescent.use_directional_spectrum` | int | `1` = directional spreading; `0` = long-crested +x sea | `1` | `0` or `1` |
| `crescent.directional_exponent` | int | Cosine-power spreading exponent | `2` | `>= 0` |
| `crescent.angle_spread` | real | Spreading width (smaller = narrower/swell) | `1.0` | in `(0, 2 pi]` |
| `crescent.jonswap_amplitude_scale` | real | Multiplier on the spectrum-derived JONSWAP Fourier amplitudes (the JONSWAP amplitude control; `wave_steepness` is ignored) | `1.0` | positive finite |
| `crescent.dt_override` | real | Explicit timestep (see Time Control); omit for T_peak/20 | automatic | positive finite |

### Derived quantities that are not inputs

`crescent.alpha`, `crescent.peak_frequency` and `crescent.peak_wavelength` are
derived from `wind_speed_10m` and `fetch`; `crescent.swell_wavenumber` has no
initializer. Setting any of them aborts the run at startup.

---

## Output Control

| Parameter | Type | Description | Default | Range |
|-----------|------|-------------|---------|-------|
| `crescent.plot_file` | string | AMReX plotfile prefix (required if plot or Tecplot output is enabled) | `""` | non-empty when enabled |
| `crescent.plot_int` | int | AMReX plotfile interval (steps); `-1` disables | `-1` | `-1` or `>= 1` |
| `crescent.tecplot_int` | int | Tecplot ASCII interval (steps); `-1` disables | `-1` | `-1` or `>= 1` |
| `crescent.verbosity` | int | Log verbosity: 0 = quiet, 1 = minimal, 2 = diagnostic | `0` | `0-2` |
| `crescent.write_omega_a_edge_diagnostics` | int | Write Omega_A y-variation edge diagnostics (L1) | `0` | `0` or `1` |
| `crescent.omega_a_edge_diagnostic_int` | int | Interval for the above (steps) | `1` | `>= 1` |

- **AMReX plotfiles** (`plot_int`): full field data, ParaView/restart-capable.
- **Tecplot** (`tecplot_int`): surface-only ASCII (X, Y, eta, u, v).
- A positive value outputs every N steps; `-1` disables; `0` is invalid.
- Runs also write CSV diagnostics (energy, including the Zakharov Hamiltonian).

---

## Checkpoint / Restart

| Parameter | Type | Description | Default | Range |
|-----------|------|-------------|---------|-------|
| `crescent.check_file` | string | Checkpoint file prefix | `"crescent_chk"` | non-empty when enabled |
| `crescent.check_int` | int | Steps between checkpoints; `-1` disables | `-1` | `-1` or `>= 1` |
| `crescent.restart` | string | Path to a checkpoint to restart from | `""` (cold start) | none |

Checkpoints now serialize the complete run configuration
(`Source/solver/CRESCENT_CheckpointConfigIO.cpp`), so a restart reproduces the
wave-initialization and nesting controls of the original run.

---

## Nested AMR (CRESCENT algorithm)

These parameters apply **only when `amr.max_level = 1`**. They are advanced; the
defaults are tuned for the standard two-level (L0 -> L1) setup. A nested L1 run
**requires** `crescent.use_crescent_algorithm = 1`, an L1 box, and a single MPI
rank.

### Region of interest & activation

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.use_crescent_algorithm` | int | 0 = single level, 1 = enable nested L1 | `0` | must be `1` for L1, `0` for L0-only |
| `crescent.L1_box_lo` | real real | L1 refined box lower (x, y) | none | required for L1 |
| `crescent.L1_box_hi` | real real | L1 refined box upper (x, y) | none | `lo < hi` |
| `crescent.roi_accurate_fraction` | real | Margin for the accurate region Omega_A | `0.25` | see below |
| `crescent.roi_spectral_fraction` | real | Margin for the spectral region Omega_S | `0.15` | `0 < spectral < accurate < 0.5` |

### Interpolation and blending mask

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.state_interpolation_order` | int | L0 -> L1 eta/phi interpolation: 1 = bilinear, 2 = quadratic | `2` | `1` or `2` |
| `crescent.derivative_interpolation_order` | int | L0 -> L1 derivative interpolation | `2` | `1` or `2` |
| `crescent.blending_steepness` | real | Has no effect; a nested run that sets it prints a warning | none | any |

The blending mask sigma is the septic tensor-product smoothstep, C^3 everywhere
including the corner seams (the former `blending_function_mode = 3`); it has no
input.

### C3 interface repair

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.c3_boundary_penalty` | real | Endpoint penalty weight lambda_b (manuscript baseline 0.01 since 2026-08-16; was 0.001) | `0.01` | non-negative finite |
| `crescent.c3_coefficient_regularization` | real | Coefficient regularization lambda_c, L2 penalty on alpha (manuscript baseline 10 since 2026-08-16; was 1) | `10.0` | non-negative finite |

PTF (mode 1) always applies the regularized C3-target interface repair in x and
y; CCF (mode 2) never uses it. The buffer widths are computed: the x and y
distances from the child box edge to Omega_S, in child cells, at least 8 cells.
A buffer that reaches into Omega_A stops the run.

### L1 coupling to L0 (PTF and CCF)

Introduced by the L1-consistency development line (synced 2026-07-17) and
qualified by the 2026-08 review campaign (synced 2026-09-04); see
`docs/CRESCENT_ALGORITHM_OVERVIEW.md` and `docs/R2_CORRECTION_FIELD_SPEC.md`.

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|------------|
| `crescent.l1_coupling_mode` | int | How the child level couples to L0: 1 = periodized total-field (PTF, production default since 2026-08-20): tendency blending on C3-periodized fields; 2 = compact correction-field (CCF): interpolated-parent background + chi-tapered compact correction confined to the child modes above the parent's retained cutoff. The historical derivative blend (0) was retired on 2026-10-02 and stops the run | `1` | `1` or `2` |

Startup constraints:
- Mode 2 (CCF) requires `state_interpolation_order == derivative_interpolation_order`.
- Nested runs require exactly one MPI rank.

---

## Example Input Files

### Minimal Stokes wave (standalone)

```bash
# amr / geometry are unprefixed in standalone runs
amr.n_cell           = 64 64 1
amr.max_level        = 0
geometry.prob_extent = 100.0 100.0 10.0
geometry.is_periodic = 1 1 0

max_step             = 100
crescent.dt          = 0.05

crescent.wave_type        = stokes
crescent.perturbation_order = 3
crescent.wave_mode_x      = 5            # k = 5 * 2 pi / 100 m: 20 m waves
crescent.wave_steepness   = 0.1

crescent.plot_file = "test"
crescent.plot_int  = 10
```

### JONSWAP wind-sea (standalone)

```bash
amr.n_cell           = 512 512 1
amr.max_level        = 0
geometry.prob_extent = 500.0 500.0 10.0
geometry.is_periodic = 1 1 0

max_step  = 4000
stop_time = 200.0
crescent.dt_override = 0.05              # JONSWAP: omit for the automatic T_peak/20

crescent.perturbation_order = 5
crescent.beta1_filter = 10.0
crescent.beta2_filter = 15.0

crescent.wave_type              = jonswap
crescent.wind_speed_10m         = 18.0
crescent.fetch                  = 15000.0
crescent.use_directional_spectrum = 1
crescent.angle_spread           = 1.2
crescent.directional_exponent   = 2

crescent.plot_file   = "steep_waves"
crescent.plot_int    = 100
crescent.tecplot_int = 20
crescent.check_file  = "crescent_chk"
crescent.check_int   = 500
```

### Nested L1 run (standalone)

```bash
amr.n_cell    = 256 256 1
amr.max_level = 1
amr.ref_ratio_vect = 2 2 1
geometry.prob_extent = 500.0 500.0 10.0
geometry.is_periodic = 1 1 0

max_step    = 1000
crescent.dt = 0.02

crescent.wave_type          = stokes
crescent.perturbation_order = 4

# Nested CRESCENT algorithm (required for L1); PTF (l1_coupling_mode = 1) is the default
crescent.use_crescent_algorithm = 1
crescent.L1_box_lo = 150.0 150.0
crescent.L1_box_hi = 350.0 350.0

crescent.plot_file = "nested"
crescent.plot_int  = 50
```

> For a coupled (ERF + CRESCENT) deck, the wave parameters are identical but
> `amr.*` / `geometry.*` are written as `crescent.amr.*` / `crescent.geometry.*`.
> See `models/ERF/Exec/DevTests/ERF_CRESCENT_Test/inputs_crescent`.

The rejection decks under `tests/inputs/` (run by `ctest` in the CRESCENT build
directory) are the executable record of the input contracts.

---

## Retired inputs

The 2026-10-02 code review removed these switches, and `ccf_band_restrict`
followed on 2026-10-03. Each one either fixed a former default or selected a
legacy or exploratory branch. A deck can still set a retired switch to the
value in the table (its supported value); any other value stops the run with a
message that names commit `f39cacf`.
`crescent.l1_coupling_mode = 0` (the derivative blend) also stops the run.

| Group | Inputs (all `crescent.`) | Accepted value |
|-------|--------------------------|----------------|
| Output | `netcdf_int`, `hdf5_int` | `-1` |
| Correction field | `ccf_band_restrict` | `1` |
| Initialization | `correct_jonswap_initialization`, `full_directional_jonswap`, `correct_stokes_wavevector`, `honor_peak_wavenumber_override` | `1` |
| Child dealiasing | `l1_square_dealias` | `1` |
| Blending mask | `blending_function_mode` | `3` |
| Blending mask | `use_sequential_blending`, `blending_correction_y_smoothing_passes` | `0` |
| Blending mask | `correct_sequential_blending` | `1` |
| Blending mask | `blending_tensor_alpha`, `blending_superellipse_gamma` | `1.0` |
| Blending mask | `blending_superellipse_p` | `4.0` |
| Blending mask | `damping_steepness` | `2.0` |
| C3 repair | `use_c3_repair`, `c3_use_interface`, `c3_apply_x`, `c3_apply_y`, `c3_directional_buffer_widths` | `1` |
| C3 repair | `c3_buffer_width` | `24` |
| C3 repair | `c3_repair_after_filter`, `c3_x_pass_y_smoothing_passes`, `c3_x_pass_corner_taper`, `c3_x_alpha_y_smoothing_passes`, `c3_x_transverse_harmonics`, `c3_y_transverse_harmonics`, `c3_x_alpha_y_projection_harmonics`, `c3_y_alpha_x_projection_harmonics` | `0` |
| C3 repair | `c3_x_alpha_y_regularization`, `c3_y_alpha_x_regularization`, `c3_x_alpha_y_mean_regularization`, `c3_y_alpha_x_mean_regularization` | `0.0` |
| Child derivatives and filters | `filter_periodic_copies_before_linear_derivatives`, `l1_filter_preserve_mean`, `enable_post_rk4_state_filter_and_blend`, `l1_correct_bilinear_interpolation`, `l1_final_dno_projection` | `1` |
| Child derivatives and filters | `linear_derivative_spectral_filter_mode`, `l1_stage_phi_z_cleanup`, `enable_post_timestep_low_pass_filter`, `enable_final_output_blended_derivative_filter`, `l1_rhs_eta_source` | `0` |
| Correction field (mode 2) | `ccf_background_tendency`, `ccf_linear_consistency`, `ccf_adjustment_consistency`, `ccf_transform_admissibility`, `ccf_w1_projection_consistency` | `1` |
| Pressure and output | `consistent_applied_pressure`, `consistent_output_diagnostics` | `1` |

---

## Parameter Validation

The solver aborts (`amrex::Abort`) on, among others:
- Missing or invalid required parameters (`crescent.dt` for `stokes` / `crescent`;
  `wind_speed_10m` / `fetch` for JONSWAP).
- `crescent.dt` given for JONSWAP, or `crescent.dt_override` given for other wave types.
- Any of `alpha`, `peak_frequency`, `peak_wavelength`, `swell_wavenumber` present.
- `perturbation_order < 1`, negative `wave_steepness`, non-positive filter widths.
- `wave_mode_x < 1` for Stokes initialization.
- Output intervals `< -1`; `plot_file` empty while plotting is enabled.
- A retired switch with a value other than its supported value (see
  [Retired inputs](#retired-inputs)), and `l1_coupling_mode = 0`.
- Out-of-range nested-AMR options (interpolation order not in {1,2},
  `l1_coupling_mode` not in {1,2}, negative or
  non-finite C3 weights, ROI fractions outside `0 < spectral < accurate < 0.5`,
  C3 buffers that reach into Omega_A, and the mode-2 constraint listed above).
- Nested L1 without `use_crescent_algorithm = 1`, without a valid L1 box, with
  `amr.max_level > 1`, or on more than one MPI rank.

---

## Appendix: Units

| Quantity | Units |
|----------|-------|
| Length (L, dx) | m |
| Time (t, dt) | s |
| Velocity (u, v, w) | m/s |
| Wavenumber (k) | rad/m |
| Frequency (omega) | rad/s |
| Wave steepness (ka) | dimensionless |
| Surface elevation (eta) | m |
| Velocity potential (phi) | m^2/s |
| Pressure (p) | Pa |

---

## Appendix: Physical Constants

Defined in `Source/shared/CRESCENT_Constants.H` (`crescent::constants`):

| Constant | Symbol | Value |
|----------|--------|-------|
| Gravity | `GRAVITY` | 9.81 m/s^2 |
| Pi | `PI` | 3.14159265... |
| 2 pi | `PI2` | 6.28318530... |
| Epsilon | `EPSILON` | 1.0e-12 |
| Small number | `SMALL_NUMBER` | 1.0e-14 |

---

## See also

- [docs/CRESCENT_ALGORITHM_OVERVIEW.md](docs/CRESCENT_ALGORITHM_OVERVIEW.md): algorithm and solver overview
- [docs/CRESCENT_CODE_WALKTHROUGH.md](docs/CRESCENT_CODE_WALKTHROUGH.md): source layout, execution paths, controls
- [README.md](README.md): quick start
- Example decks under `Exec/` and `test_case/`

*Questions: the CRESCENT team at LLNL.*

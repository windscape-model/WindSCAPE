# CRESCENT Input Parameter Reference

Complete reference for the input parameters read by the CRESCENT wave solver
(High-Order Spectral, HOS). Generated from the parameter readers in
`Source/solver/CRESCENT_SetupRuntimeInputs.cpp`,
`Source/solver/CRESCENT_SetupParameters.cpp`, and
`Source/nested/CRESCENT_NestedAlgorithmInputs.cpp`; defaults are taken from the
state structs (`CRESCENT_RuntimeControl.H`, `CRESCENT_WaveState.H`,
`CRESCENT_NestedConfig.H`).

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

---

## Time Control

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `max_step` | int | Maximum number of timesteps (global, unprefixed) | `INT_MAX` (∞) |
| `stop_time` | real | Simulation end time, s (global, unprefixed) | `REAL_MAX` (∞) |
| `start_time` | real | Simulation start time, s (global, unprefixed) | `0.0` |
| `crescent.dt` | real | Fixed timestep size, s | none (must be set) |
| `crescent.dt_override` | real | Advanced override of `dt`; `-1.0` = use `crescent.dt` | `-1.0` |

The run ends at whichever of `max_step` / `stop_time` is reached first.

---

## Domain (AMR & Geometry)

Standard AMReX `AmrCore` parameters. In **standalone** runs these are unprefixed;
in **coupled** runs prepend `crescent.` (e.g. `crescent.amr.n_cell`).

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `amr.n_cell` | int int int | Grid cells (x, y, z) | Required |
| `amr.max_level` | int | Maximum AMR level (0 = single level) | `0` |
| `amr.max_grid_size` | int | Maximum grid (box) size | AMReX default |
| `amr.blocking_factor` | int | Grid blocking factor | AMReX default |
| `geometry.prob_lo` | real real real | Lower domain corner, m | `0 0 0` |
| `geometry.prob_hi` | real real real | Upper domain corner, m | — |
| `geometry.prob_extent` | real real real | Domain size, m (alternative to `prob_hi`) | — |
| `geometry.is_periodic` | int int int | Periodic BC flags (x, y, z) | Required (use `1 1 0`) |

**Notes**
- The HOS spectral method is single-level on L0; set `amr.max_level = 0` unless
  you are using the nested CRESCENT algorithm (see below).
- Use periodic in x and y for ocean waves: `geometry.is_periodic = 1 1 0`.
- `n_cell` in x/y should be powers of 2 for FFT efficiency.
- Domain size sets the wavenumber resolution Δk = 2π/L (computed automatically).

---

## Spectral / HOS Solver

| Parameter | Type | Description | Default | Range |
|-----------|------|-------------|---------|-------|
| `crescent.perturbation_order` | int | HOS expansion order M | `3` | `≥ 1` |
| `crescent.beta1_filter` | real | Spectral filter parameter 1 | `8.0` | — |
| `crescent.beta2_filter` | real | Spectral filter parameter 2 | `12.0` | — |
| `crescent.peak_wavenumber_override` | real | Peak wavenumber for the filter; `-1.0` = compute automatically | `-1.0` | — |

### Perturbation Order

`perturbation_order` (M) trades accuracy for cost. Higher M captures steeper,
more nonlinear waves:

| M | Regime | Typical use |
|---|--------|-------------|
| 1 | Linear | Testing only |
| 2 | Weakly nonlinear | Mild seas |
| 3 | Nonlinear | **Recommended default** |
| 4–5 | Highly nonlinear | Steep / near-breaking waves |
| 6+ | Extreme | Research |

> Wavenumber spacing (Δkₓ, Δk_y) is **computed internally** from the domain size;
> it is not an input parameter.

---

## Wave Initialization

The initialization model is chosen with `crescent.wave_type`.

### Common

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `crescent.wave_type` | string | `"jonswap"`, `"stokes"`, `"monochromatic"`, or `"focused_group"` | `"jonswap"` |
| `crescent.wave_steepness` | real | Wave steepness ka | `0.1` |
| `crescent.initial_angle` | real | Initial/mean wave direction, rad | `0.0` |
| `crescent.wave_angle` | real | Stokes propagation angle, rad (0 = +x, π/2 = +y) | `0.0` |

### Stokes waves (`crescent.wave_type = stokes`)

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|-----------|
| `crescent.swell_wavenumber` | real | Primary wavenumber k, rad/m | `1.0` | — |
| `crescent.wave_mode_x` | int | Wave mode index (sets wavenumber magnitude) | `5` | `≥ 1` |

Wavelength relation: λ = 2π / k. Amplitude from steepness: a = ka / k.

### JONSWAP spectrum (`crescent.wave_type = jonswap`)

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|-----------|
| `crescent.wind_speed_10m` | real | Wind speed at 10 m, m/s | `0.0` | `> 0` required |
| `crescent.fetch` | real | Fetch length, m | `0.0` | `> 0` required |
| `crescent.jonswap_random_seed` | long | Phase seed for reproducible init; `-1` = AMReX default seeding | `12345` | `≥ -1` |
| `crescent.use_directional_spectrum` | int | Enable directional spreading | `1` | `0` or `1` |
| `crescent.directional_exponent` | int | Cosine spreading exponent | `2` | — |
| `crescent.angle_spread` | real | Spreading width (smaller = narrower/swell) | `1.0` | — |

### Advanced spectrum overrides

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `crescent.alpha` | real | Phillips/JONSWAP α (0 = derive from wind & fetch) | `0.0` |
| `crescent.peak_frequency` | real | Spectral peak frequency, rad/s | `1.0` |
| `crescent.peak_wavelength` | real | Spectral peak wavelength, m | `0.0` |

---

## Output Control

| Parameter | Type | Description | Default | Range |
|-----------|------|-------------|---------|-------|
| `crescent.plot_file` | string | AMReX plotfile prefix (required if plotting enabled) | `""` | — |
| `crescent.plot_int` | int | AMReX plotfile interval (steps); `-1` disables | `-1` | `-1` or `≥ 1` |
| `crescent.tecplot_int` | int | Tecplot ASCII interval (steps); `-1` disables | `-1` | `-1` or `≥ 1` |
| `crescent.netcdf_int` | int | NetCDF interval, **not implemented** (`> 0` aborts) | `-1` | `-1` |
| `crescent.hdf5_int` | int | HDF5 interval, **not implemented** (`> 0` aborts) | `-1` | `-1` |
| `crescent.verbosity` | int | Log verbosity: 0=quiet, 1=minimal, 2=diagnostic | `0` | `0–2` |
| `crescent.write_omega_a_edge_diagnostics` | int | Write Ω_A y-variation edge diagnostics (L1) | `0` | `0` or `1` |
| `crescent.omega_a_edge_diagnostic_int` | int | Interval for the above (steps) | `1` | `≥ 1` |

- **AMReX plotfiles** (`plot_int`): full field data, ParaView/restart-capable.
- **Tecplot** (`tecplot_int`): surface-only ASCII (X, Y, η, u, v).
- A positive value outputs every N steps; `-1` disables; `0` is invalid.

---

## Checkpoint / Restart

| Parameter | Type | Description | Default | Range |
|-----------|------|-------------|---------|-------|
| `crescent.check_file` | string | Checkpoint file prefix | `"crescent_chk"` | — |
| `crescent.check_int` | int | Steps between checkpoints; `-1` disables | `-1` | `-1` or `≥ 1` |
| `crescent.restart` | string | Path to a checkpoint to restart from | `""` (cold start) | — |

---

## Nested AMR (CRESCENT algorithm)

These parameters apply **only when `amr.max_level ≥ 1`**. They are advanced; the
defaults are tuned for the standard two-level (L0→L1) setup. A nested L1 run
**requires** `crescent.use_crescent_algorithm = 1` and an L1 box.

### Region of interest & activation

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|-----------|
| `crescent.use_crescent_algorithm` | int | 0 = single level, 1 = enable nested L1 | `0` | must be `1` for L1 |
| `crescent.L1_box_lo` | real real | L1 refined box lower (x, y) | — | required for L1 |
| `crescent.L1_box_hi` | real real | L1 refined box upper (x, y) | — | `lo < hi` |
| `crescent.roi_accurate_fraction` | real | Margin for the accurate region Ω_A | `0.25` | — |
| `crescent.roi_spectral_fraction` | real | Margin for the spectral region Ω_S | `0.15` | — |

### Interpolation & blending

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|-----------|
| `crescent.state_interpolation_order` | int | L0→L1 η/φ interp: 1=bilinear, 2=quadratic | `2` | `1` or `2` |
| `crescent.derivative_interpolation_order` | int | L0→L1 derivative interp | `2` | `1` or `2` |
| `crescent.blending_function_mode` | int | 0=Euclidean, 1=sequential, 2=harmonic, 3=tensor, 4=harmonic+septic, 5=superellipse | `0` | `0–5` |
| `crescent.use_sequential_blending` | int | 0=2D Euclidean, 1=sequential 1D | `0` | `0` or `1` |
| `crescent.blending_tensor_alpha` | real | Tensor-product σ shape (mode 3) | `1.0` | `> 0` |
| `crescent.blending_superellipse_gamma` | real | Superellipse side-ramp exponent (mode 5) | `1.0` | `> 0` |
| `crescent.blending_superellipse_p` | real | Superellipse exponent (mode 5) | `4.0` | `> 0` |
| `crescent.blending_correction_y_smoothing_passes` | int | 1-2-1 passes on blend correction along y | `0` | `≥ 0` |
| `crescent.damping_steepness` | real | Damping-path shape; ignored when C3 repair active | `2.0` | — |
| `crescent.blending_steepness` | real | Accepted but currently unused by active σ builders | `2.0` | — |

### C3 interface repair

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|-----------|
| `crescent.use_c3_repair` | int | 0 = damping path, 1 = C3 repair | `1` | — |
| `crescent.c3_apply_x` | int | Apply x-pass repair | `1` | `0` or `1` |
| `crescent.c3_apply_y` | int | Apply y-pass repair | `1` | `0` or `1` |
| `crescent.c3_use_interface` | int | 0=cell-center, 1=interface-based | `1` | — |
| `crescent.c3_buffer_width` | int | Repair buffer width (grid points) | `24` | — |
| `crescent.c3_boundary_penalty` | real | Boundary magnitude constraint weight | `0.001` | — |
| `crescent.c3_coefficient_regularization` | real | L2 penalty on repair coefficients | `1.0` | — |
| `crescent.c3_repair_after_filter` | int | Second C3 repair after L1 filtering | `0` | `0` or `1` |
| `crescent.c3_x_pass_y_smoothing_passes` | int | 1-2-1 passes on x-pass correction along y | `0` | `≥ 0` |
| `crescent.c3_x_pass_corner_taper` | int | Taper x-pass in y-buffers | `0` | `0` or `1` |
| `crescent.c3_x_alpha_y_smoothing_passes` | int | 1-2-1 passes on x-pass coefficients along y | `0` | `≥ 0` |
| `crescent.c3_x_alpha_y_regularization` | real | x-pass coefficient smoothness penalty (y) | `0.0` | `≥ 0` |
| `crescent.c3_y_alpha_x_regularization` | real | y-pass coefficient smoothness penalty (x) | `0.0` | `≥ 0` |
| `crescent.c3_x_alpha_y_mean_regularization` | real | x-pass transverse-mean deviation penalty | `0.0` | `≥ 0` |
| `crescent.c3_y_alpha_x_mean_regularization` | real | y-pass transverse-mean deviation penalty | `0.0` | `≥ 0` |
| `crescent.c3_x_transverse_harmonics` | int | Fourier harmonics for x-pass coeffs (0=linewise) | `0` | `≥ 0` |
| `crescent.c3_y_transverse_harmonics` | int | Fourier harmonics for y-pass coeffs (0=linewise) | `0` | `≥ 0` |
| `crescent.c3_x_alpha_y_projection_harmonics` | int | Fourier projection of x-pass coeffs along y | `0` | `≥ 0` |
| `crescent.c3_y_alpha_x_projection_harmonics` | int | Fourier projection of y-pass coeffs along x | `0` | `≥ 0` |

> When `use_c3_repair = 1`, `c3_apply_x` and `c3_apply_y` cannot both be 0.

### Linear-derivative filtering & diagnostics

| Parameter | Type | Description | Default | Constraint |
|-----------|------|-------------|---------|-----------|
| `crescent.filter_periodic_copies_before_linear_derivatives` | int | 0=post-C3 copies, 1=filtered copies | `1` | `0` or `1` |
| `crescent.linear_derivative_spectral_filter_mode` | int | 0=2/3 cutoff, 1=+Hou-Li, 2=+Hou-Li+low-pass | `0` | `0–2` |
| `crescent.enable_post_timestep_low_pass_filter` | int | Low-pass after surface update | `0` | `0` or `1` |
| `crescent.enable_post_rk4_state_filter_and_blend` | int | Run nested post-RK4 filter/blend | `1` | `0` or `1` |
| `crescent.enable_final_output_blended_derivative_filter` | int | Post-blend FFT filtering of output derivatives | `0` | `0` or `1` |
| `crescent.enable_initial_c3_effect_diagnostic` | int | First-step C3 effect diagnostic | `1` | `0` or `1` |
| `crescent.enable_nested_initialization_validation` | int | Nested initialization validation bundle | `1` | `0` or `1` |

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
crescent.swell_wavenumber = 0.314159     # 20 m waves
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
crescent.dt = 0.05

crescent.perturbation_order = 5
crescent.beta1_filter = 10.0
crescent.beta2_filter = 15.0

crescent.wave_type              = jonswap
crescent.wind_speed_10m         = 18.0
crescent.fetch                  = 15000.0
crescent.wave_steepness         = 0.22
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
geometry.prob_extent = 500.0 500.0 10.0
geometry.is_periodic = 1 1 0

max_step    = 1000
crescent.dt = 0.02

crescent.wave_type          = stokes
crescent.perturbation_order = 4

# Nested CRESCENT algorithm (required for L1)
crescent.use_crescent_algorithm = 1
crescent.L1_box_lo = 150.0 150.0
crescent.L1_box_hi = 350.0 350.0
crescent.use_c3_repair = 1

crescent.plot_file = "nested"
crescent.plot_int  = 50
```

> For a coupled (ERF + CRESCENT) deck, the wave parameters are identical but
> `amr.*` / `geometry.*` are written as `crescent.amr.*` / `crescent.geometry.*`.
> See `models/ERF/Exec/DevTests/ERF_CRESCENT_Test/inputs_crescent`.

---

## Parameter Validation

The solver aborts (`amrex::Abort`) on, among others:
- Missing required parameters (`crescent.dt`; `wind_speed_10m`/`fetch` for JONSWAP).
- `perturbation_order < 1`, `wave_steepness < 0`.
- `wave_mode_x < 1` for Stokes initialization.
- Output intervals `< -1`; `plot_file` empty while plotting is enabled.
- `netcdf_int > 0` or `hdf5_int > 0` (not implemented).
- Out-of-range nested-AMR options (interpolation order ∉ {1,2},
  `blending_function_mode` ∉ {0..5}, negative regularizations, etc.).
- Nested L1 without `use_crescent_algorithm = 1` or without a valid L1 box.

---

## Appendix: Units

| Quantity | Units |
|----------|-------|
| Length (L, dx) | m |
| Time (t, dt) | s |
| Velocity (u, v, w) | m/s |
| Wavenumber (k) | rad/m |
| Frequency (ω) | rad/s |
| Wave steepness (ka) | dimensionless |
| Surface elevation (η) | m |
| Velocity potential (φ) | m²/s |
| Pressure (p) | Pa |

---

## Appendix: Physical Constants

Defined in `Source/shared/CRESCENT_Constants.H` (`crescent::constants`):

| Constant | Symbol | Value |
|----------|--------|-------|
| Gravity | `GRAVITY` | 9.81 m/s² |
| Pi | `PI` | 3.14159265… |
| 2π | `PI2` | 6.28318530… |
| Epsilon | `EPSILON` | 1.0e-12 |
| Small number | `SMALL_NUMBER` | 1.0e-14 |

---

## See also

- [docs/CRESCENT_ALGORITHM_OVERVIEW.md](docs/CRESCENT_ALGORITHM_OVERVIEW.md): algorithm and solver overview
- [README.md](README.md): quick start
- Example decks under `Exec/` and `test_case/`

*Questions: the CRESCENT team at LLNL.*

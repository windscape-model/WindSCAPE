# Wind over waves

These decks simulate air over a prescribed monochromatic water wave on a moving terrain-following
mesh, with the split-explicit (acoustic sub-step) time integration of ERF. They are verification
cases of the manuscript "A split-explicit moving-mesh method for wind over ocean waves in the
Energy Research and Forecasting (ERF) model". Each case runs on a workstation.

## Build

From `Exec/`:

    make -j8 USE_NETCDF=FALSE                             # all cases except the Stokes layer
    make -j8 USE_NETCDF=FALSE USE_TERRAIN_VELOCITY=TRUE   # Stokes layer: executable ERF3d.*.TV.ex

The cases do not need NetCDF. `USE_TERRAIN_VELOCITY = TRUE` (CMake:
`-DERF_ENABLE_TERRAIN_VELOCITY=ON`) sets the velocity of a no-slip lower wall to the orbital
velocity of the wave. Without it, a no-slip lower wall is at rest.

## Cases

| Deck | Case |
|---|---|
| `inputs_inviscid_wave` | traveling wave under a slip lid, no viscosity |
| `inputs_reflection_free_wave` | standing wave started from rest in a deep domain, no off-centering and no divergence damping |
| `inputs_stratified_wave` | the reflection-free wave over a stably stratified base state |
| `inputs_stokes_layer` | viscous Stokes layer on the no-slip wave surface (TV build) |

Run `inputs_stratified_wave` from a folder that contains its sounding file `input_sounding_stratified`.
For a study of the time step and the acoustic sub-step, change `erf.fixed_dt` and
`erf.fixed_mri_dt_ratio` (M); the sub-step is `erf.fixed_dt / M`.

## Inputs of the moving wave

| Input | Default | Meaning |
|---|---|---|
| `prob.custom_terrain_type` | | `MovingSineWave`: h = a sin(kx - omega t); `MovingSineWaveRest`: h = a sin(kx) (1 - cos(omega t)); `MovingSineWaveRamp`: h = a R(t) sin(kx - omega t), with R(t) = (1 - cos(omega t))/2 for omega t < pi and 1 after |
| `prob.Ampl` | 0 | amplitude a [m] |
| `prob.wavelength` | 100 | wavelength [m], k = 2 pi / wavelength |
| `prob.wave_c` | 0 | phase speed c [m/s]; omega = c k if c > 0, else sqrt(g k) |
| `prob.init_wave_field` | 1 | 1: initial velocity and pressure of the potential flow over the wave; 0: air at rest |

The time integration on the moving mesh uses these evaluation times in the acoustic sub-steps. The
defaults are the scheme of the manuscript; the other values only show its conditions.

| Input | Default | Other value |
|---|---|---|
| `erf.substep_metric_sampling` | `midpoint` (moving mesh): metric terms of the pressure gradient and operator Jacobians at the sub-step midpoint | `endpoint`: at the sub-step start |
| `erf.grid_velocity_sampling` | `DiscreteDisplacement`: displacement over the sub-step | `IntervalStart`: rate over the preceding interval |
| `erf.slow_geometry_time` | `StageTime`: slow tendency on the mesh at the stage time | `StepStart`: on the mesh at the step start |
| `erf.substep_volume_sampling` | `Endpoint`: cell volume at the sub-step start on the sub-step-start state | `Midpoint`: mean of the start and end volumes |
| `erf.beta_d` | -1: 0.1, or 0 with midpoint metric terms and `erf.beta_s = 0` | divergence damping coefficient |

# CRESCENT: A block-structured, adaptive-mesh framework for ocean surface wave modeling

CRESCENT is a block-structured, adaptive-mesh research software framework for
simulating ocean surface waves, built on the AMReX adaptive mesh refinement
library. It models the time evolution of nonlinear wave fields and supports
phase-resolved numerical studies of wave dynamics, including wave groups and
broadband sea states.

The current release provides a high-order spectral (HOS) free-surface wave
solver with two methods for one-way nested refinement of a non-periodic child
region inside a periodic parent: the periodized total-field method (PTF), which
is the default, and the compact correction-field method (CCF). The
executable target is `crescent`, and inputs use the `crescent.` parameter
prefix.

## Capabilities

- Single-domain core solver for nonlinear free-surface waves represented by
  surface elevation `eta` and surface potential `phi`, with a selectable HOS
  order (`crescent.perturbation_order`).
- Nested L1 solver path with one-way L0-to-L1 coupling for non-periodic
  refinement regions. `crescent.l1_coupling_mode` selects the coupling method:
  - `1`, the periodized total-field method (PTF, default), evolves the total
    child state. It makes periodic scratch copies of the child fields for the
    spectral operators and blends the child and parent tendencies through a
    smooth mask.
  - `2`, the compact correction-field method (CCF), evolves a tapered
    correction to the interpolated parent field, only at the wavenumbers
    above the parent's cutoff. The correction is periodic by
    construction, so this method needs no periodization.
- Integer refinement ratios of two or more between parent and child, per-level
  spectral filter control (`crescent.beta1_filter_l1`), and an optional
  Dommermuth-style nonlinear ramp (`crescent.adjustment_periods`) for the
  initialization of steep seas (single-domain core and both nested methods).
- AMReX-based mesh, distribution, plotfile, checkpoint, MPI, and GPU-enabled
  infrastructure.
- Wave initializations: linear Stokes, JONSWAP broadband seas with fixed
  random seeds, and the fifth-order crescent (class-II horseshoe) case from
  which the solver takes its name (`crescent.wave_type = "crescent"`).
- Curated sample input cases for single-domain and nested wave simulations.
- Public non-owning exchange views of the surface state for external solvers.

## Build

Prerequisites:

- CMake 3.14 or newer
- C++17 compiler
- MPI C++ compiler wrapper
- FFTW3 for CPU builds (the AMReX FFT backend; GPU builds use the vendor FFT
  library, for example cuFFT)
- AMReX through the included submodule (pinned: release 25.01, commit
  `92d35c2c8`, 2025-01-18), or an externally supplied AMReX build

Typical build from the repository root:

```bash
git submodule update --init --recursive
cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build --target crescent -j4
```

The included AMReX submodule is enabled by default through
`CRESCENT_USE_INTERNAL_AMREX=ON`. To link against an external AMReX installation,
configure with `-DCRESCENT_USE_INTERNAL_AMREX=OFF` and provide the appropriate
AMReX package path through `AMReX_DIR` or `CMAKE_PREFIX_PATH` for the target
platform.

### GPU build

Both nested methods run on GPUs. For CUDA:

```bash
cmake -S . -B build_cuda \
      -DCMAKE_BUILD_TYPE=Release \
      -DCRESCENT_ENABLE_CUDA=ON \
      -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build_cuda --target crescent -j4
```

`CRESCENT_ENABLE_HIP` and `CRESCENT_ENABLE_SYCL` select the other AMReX
backends; CUDA is the tested GPU path. On a GPU, all per-cell work runs as
device kernels, transforms use `amrex::FFT::R2C` (cuFFT on CUDA), the PTF
interface repair runs as a batched per-grid-line device kernel, and output
stages through pinned host memory. GPU and CPU runs of the same deck agree to
round-off, not bit for bit.

## Run

Curated sample inputs live in [test_case](test_case/):

```bash
ROOT=$(pwd)
CASE="$ROOT/runs/sample_stokes_$(date +%Y-%m-%d_%H%M%S)"
mkdir -p "$CASE"
cp test_case/inputs_stokes_test "$CASE/inputs"
(cd "$CASE" && "$ROOT/build/Exec/crescent" inputs > run.log 2>&1)
```

Replace the copied deck with `inputs_CRESCENT_L0_only` or
`inputs_nested_jonswap_init` for the other samples. Single-domain MPI
execution uses the same fresh-directory pattern, with the invocation replaced
by:

```bash
mpiexec -n 4 "$ROOT/build/Exec/crescent" inputs > run.log 2>&1
```

Input parameters use the `crescent.` prefix.

Nested runs currently use a fixed AMR hierarchy after initialization because
the spectral and nested scratch workspaces are built for the initialized level
layout. Set `amr.regrid_int` larger than the run length for nested cases.
Nested runs are single-rank; single-level runs support MPI.

## Nested Refinement

A nested run refines a rectangular child window (level L1) inside the periodic
parent (level L0) and couples them one way: L0 advances independently, and L1
reads L0 state and tendencies without feeding back. The child is enabled and
configured with:

```ini
crescent.use_crescent_algorithm = 1      # enable the nested L1 path
crescent.l1_coupling_mode = 1            # 1 = PTF (default), 2 = CCF
crescent.L1_box_lo = 80.0 80.0           # child window, parent coordinates
crescent.L1_box_hi = 144.0 144.0
amr.max_level = 1
amr.ref_ratio_vect = 2 2 1               # integer ratio >= 2 (equal in x and y)
```

PTF is the default: the child evolves the total field and can depart from the
parent. CCF keeps the parent's resolved waves and adds the shorter waves above
the parent's cutoff: its correction holds only the child modes above the
parent's retained cutoff, and the child's shared band follows the
interpolated parent. The algorithm and the two methods are described in
[docs/CRESCENT_ALGORITHM_OVERVIEW.md](docs/CRESCENT_ALGORITHM_OVERVIEW.md).

Minimal single-domain example:

```ini
max_step = 100
stop_time = 10.0

amr.max_level = 0
amr.n_cell = 128 128 1
geometry.prob_extent = 200.0 200.0 10.0
geometry.is_periodic = 1 1 0

crescent.dt = 0.02
crescent.perturbation_order = 3
crescent.wave_type = "stokes"
crescent.wave_mode_x = 5
crescent.wave_steepness = 0.1

crescent.plot_file = "waves"
crescent.plot_int = 20
crescent.tecplot_int = 20
```

For JONSWAP seas (`crescent.wave_type = "jonswap"`), `crescent.wind_speed_10m`
and `crescent.fetch` set the spectrum, and `crescent.jonswap_random_seed` sets
the random phases. `crescent.use_directional_spectrum = 0` selects a
long-crested spectrum that travels in `+x`, and
`crescent.jonswap_amplitude_scale` multiplies all Fourier amplitudes.
`crescent.wave_steepness` applies only to the Stokes and crescent
initializers. A positive `crescent.peak_wavenumber_override` replaces the
initializer's peak wavenumber in the spectral filters.

Some options of earlier versions are retired (for example the derivative-blend
coupling, the alternative blending masks, the pre-2026-08 JONSWAP and Stokes
constructors, and `crescent.ccf_band_restrict`). An input deck can still give a
retired option its supported value; any other value stops the run with a
message that names commit `f39cacf` of the development repository, where those
options remain.

## Output

CRESCENT can write AMReX plotfiles and Tecplot-style surface files. Output is
controlled by:

- `crescent.plot_file`: output path prefix.
- `crescent.plot_int`: AMReX plotfile interval, with `-1` disabling plotfiles.
- `crescent.tecplot_int`: Tecplot surface output interval, with `-1` disabling
  Tecplot output.
- `crescent.write_omega_a_edge_diagnostics`: optional nested-region diagnostic
  CSV output.

Nested runs write their energy history to `crescent_energy.csv` in the working
directory. With `crescent.write_omega_a_edge_diagnostics = 1`, they also write
the per-step diagnostic files `omega_a_edge_diagnostic.csv` and
`nested_oscillation_diagnostic.csv` to the directory of `crescent.plot_file`.

## Checkpoint And Restart

Checkpointing is controlled with AMReX-style runtime inputs:

```ini
crescent.check_file = "crescent_chk"
crescent.check_int = 50
```

This writes checkpoint directories such as `crescent_chk00050`. To restart from a
checkpoint, run with a compatible input file that includes:

```ini
crescent.restart = "crescent_chk00050"
```

Restart restores the free-surface state, the step and time counters, the level
geometry and layout, and the wave, spectral, and nested evolution settings.
When the restart deck gives other values for these settings, the values in the
checkpoint apply; the run length and the output controls come from the restart
deck. Restart supports single-domain and nested runs. The solver reads only
checkpoints of its own format (format 10).

## External Solver Exchange

CRESCENT exposes a lightweight interface for external solvers through
`crescent::external_exchange::SurfaceExchangeProvider`. The interface returns
non-owning AMReX views of the current free-surface state for a requested level:

- `eta`: surface elevation.
- `phi`: surface potential.
- `pressure`: pressure field.
- `geometry`, `box_array`, and `distribution_mapping`: AMReX metadata needed by
  external consumers.

This API is infrastructure only. It does not impose a coupling algorithm or
modify the core wave solver.

## Documentation

- [docs/CRESCENT_ALGORITHM_OVERVIEW.md](docs/CRESCENT_ALGORITHM_OVERVIEW.md): public algorithm and architecture overview, including both nested coupling methods.
- [docs/CRESCENT_CODE_WALKTHROUGH.md](docs/CRESCENT_CODE_WALKTHROUGH.md): source-level execution and ownership guide.
- [docs/R2_CORRECTION_FIELD_SPEC.md](docs/R2_CORRECTION_FIELD_SPEC.md): equation-level CCF specification.
- [docs/README.md](docs/README.md): release documentation index.
- [test_case/README.md](test_case/README.md): curated sample cases, with run instructions.

Internal planning notes and development logs are not part of the release
documentation set.

## Contributors And Provenance

Contributor attribution and development-assistance provenance are recorded in
[CONTRIBUTORS.md](CONTRIBUTORS.md).

## License

CRESCENT is distributed under the BSD 3-Clause License; see [LICENSE](LICENSE).
It was developed at Lawrence Livermore National Laboratory. Version 1.0 is
released under LLNL release number LLNL-CODE-2024694; later versions inherit
this license.

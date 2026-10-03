# CRESCENT: A block-structured, adaptive-mesh framework for ocean surface wave modeling

CRESCENT is a block-structured, adaptive-mesh research software framework for
simulating ocean surface waves, built on the AMReX adaptive mesh refinement
library. It models the time evolution of nonlinear wave fields and supports
phase-resolved numerical studies of wave dynamics, including wave groups and
broadband sea states.

The current release provides a high-order spectral free-surface wave solver,
with a nested-grid Repair-Compute-Blend path for localized high-resolution wave
dynamics. The executable target is `crescent`, and inputs use the `crescent.`
parameter prefix.

## Capabilities

- Single-domain core solver for nonlinear free-surface waves represented by
  surface elevation `eta` and surface potential `phi`.
- Nested L1 solver path with one-way L0-to-L1 coupling, C3 periodization, and
  sigma blending for non-periodic refinement regions.
- AMReX-based mesh, distribution, plotfile, checkpoint, MPI, and GPU-enabled
  infrastructure.
- Curated sample input cases for single-domain and nested wave simulations.
- Public non-owning exchange views for future coupling to external solvers.

## Build

Prerequisites:

- CMake 3.14 or newer
- C++17 compiler
- MPI C++ compiler wrapper
- AMReX through the included submodule, or an externally supplied AMReX build
- FFT support through AMReX configuration

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

## Run

Curated sample inputs live in [test_case](test_case/):

```bash
./build/Exec/crescent test_case/inputs_stokes_test
./build/Exec/crescent test_case/inputs_CRESCENT_L0_only
./build/Exec/crescent test_case/inputs_nested_jonswap_init
```

MPI execution uses the same input files:

```bash
mpiexec -n 4 ./build/Exec/crescent test_case/inputs_stokes_test
```

Input parameters use the `crescent.` prefix.

Nested runs currently use a fixed AMR hierarchy after initialization because
the spectral and C3 workspaces are built for the initialized level layout. Set
`amr.regrid_int` larger than the run length for nested cases.

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

## Output

CRESCENT can write AMReX plotfiles and Tecplot-style surface files. Output is
controlled by:

- `crescent.plot_file`: output path prefix.
- `crescent.plot_int`: AMReX plotfile interval, with `-1` disabling plotfiles.
- `crescent.tecplot_int`: Tecplot surface output interval, with `-1` disabling
  Tecplot output.
- `crescent.write_omega_a_edge_diagnostics`: optional nested-region diagnostic
  CSV output.

NetCDF and HDF5 interval keys are reserved, but those writers are not currently
implemented.

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

Restart preserves the free-surface state, step/time counters, level geometry and
layout metadata, runtime reporting state, nested configuration, and enough
metadata to rebuild derived fields and scratch workspaces safely. The restart
path supports both single-domain and nested configurations.

## External Solver Exchange

CRESCENT exposes a lightweight future-coupling interface through
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

- [docs/CRESCENT_ALGORITHM_OVERVIEW.md](docs/CRESCENT_ALGORITHM_OVERVIEW.md): public algorithm and architecture overview.
- [docs/README.md](docs/README.md): release documentation index.
- [test_case/README.md](test_case/README.md): curated sample cases.

Internal planning notes and development logs are not part of the release
documentation set.

## Contributors And Provenance

Contributor attribution and development-assistance provenance are recorded in
[CONTRIBUTORS.md](CONTRIBUTORS.md).

## License

License and release metadata will be added through the LLNL software release
approval process before public distribution.

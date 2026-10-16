# WINDSCAPE: Wind over a SCAPE, a Single-executable Coupled AMReX Peer Environment

> **Project:** Coupled atmosphere–wave modelling for offshore wind and air–sea interaction
> **Funding:** LLNL Laboratory Directed Research and Development (LDRD), project 25-LW-087
> **Institution:** Lawrence Livermore National Laboratory
> **Contact:** Tianyi Li (tianyili1993@gmail.com)

WINDSCAPE is an open, two-way **coupled atmosphere–wave model** for the phase-resolved marine
atmospheric boundary layer, built on the [AMReX](https://amrex-codes.github.io/amrex/)
block-structured adaptive-mesh framework. It couples:

- **ERF** (Energy Research and Forecasting): a compressible large-eddy / mesoscale atmospheric solver, and
- **CRESCENT**: a nested, phase-resolved High-Order Spectral (HOS) nonlinear free-surface wave solver,

**in process, inside a single executable.** A top-level *Coupler* builds ERF and CRESCENT as
equal-weight peers that share one build of AMReX, and installs the coupling exchange into ERF's
coarse time loop as a per-step hook (`ERF::Evolve` hosts the loop). The coupling is
two-way and tight: every atmospheric step, ERF's surface pressure forces the waves, and the
evolving wave surface drives ERF's terrain-following lower boundary through a moving fitted mesh.

The target application is **offshore wind** and air–sea interaction, where the wind in the lowest
tens of metres is strongly modulated by a moving, phase-resolved sea surface that bulk-roughness
parameterisations and spectral wave models cannot represent.

---

## Repository layout

```
WINDSCAPE/
  CMakeLists.txt          top-level superbuild (CMake ExternalProject)
  Coupler/                the coupling driver (WINDSCAPE's own code)
    main.cpp              amrex::Initialize -> build ERF + CRESCENT -> run -> Finalize
    Coupler.{H,cpp}       ERF <-> CRESCENT transfer kernels
    CMakeLists.txt
  models/ERF/             vendored ERF (coupling-agnostic; thin coupling hooks, see UPSTREAM.md)
  models/CRESCENT/        vendored CRESCENT (full solver, including nesting; coupling adapter)
  shared/amrex/           one vendored AMReX, shared by both models (two local patches)
  tests/                  standalone unit tests (e.g. per-instance geometry)
  Documentation/          WindSCAPE docs (CI design and test tiers: CI.md)
  .github/workflows/      active CI (GitHub Actions: superbuild + smoke suite)
  .gitlab-ci.yml, .gitlab/  dormant LLNL LC GitLab pipeline (see Documentation/CI.md)
```

This is a **self-contained superproject**: ERF, CRESCENT, and AMReX are vendored as plain files
(no submodules), so a checkout is fully buildable offline and the provenance of every component
is recorded (see [Components and versions](#components-and-versions)).

## Building

Requirements: CMake ≥ 3.20, an MPI installation, FFTW3 (single- and double-precision, used by
CRESCENT's spectral solver through AMReX's FFT backend), and a C++17 compiler (GPU builds are
optional via AMReX's CUDA/HIP/SYCL backends).

The top-level build is a CMake `ExternalProject` superbuild that builds and installs the shared
AMReX once, then builds ERF and CRESCENT as libraries against it, then links the coupled executable:

```bash
cmake -S . -B build/superbuild
cmake --build build/superbuild -j
# -> build/superbuild/coupler-build/windscape
```

(`build/` is a generated artifact directory and is git-ignored.)

## Running a coupled case

A reference coupled case lives under the ERF problem directory. Run the coupled executable on the
ERF input deck, which selects the wave-driven moving terrain and points to the wave deck:

```bash
cd models/ERF/Exec/DevTests/ERF_CRESCENT_Test
mpiexec -np N /path/to/build/superbuild/coupler-build/windscape inputs_erf
```

Key coupling controls in `inputs_erf` (`erf.crescent_orbital_bc` is shown with its default; the
deck does not set it):

```
erf.terrain_type         = MovingFittedMesh   # wave-driven, terrain-following lower boundary
prob.custom_terrain_type = "CRESCENTWaves"
erf.crescent_coupling      = true
erf.crescent_input_file    = "inputs_crescent"
erf.crescent_dt_ratio      = 5                  # CRESCENT sub-steps per ERF coarse step
erf.crescent_grid_ratio    = 5                  # ERF:CRESCENT horizontal grid ratio
erf.crescent_orbital_bc    = true               # default: push the wave orbital u for the moving-wall BC
```

With a no-slip lower wall (`zlo.type = "NoSlipWall"`) and molecular viscosity, the coupled run
also imposes CRESCENT's orbital x-velocity at the surface through ERF's gated moving-wall BC
(`ERF_ENABLE_TERRAIN_VELOCITY`, on in the superbuild): see the `inputs_erf_orbital` variant
deck. With a slip wall the BC is inert and results are unchanged.

The wave solver reads its own parameters under the `crescent.` prefix from `inputs_crescent`. A nested
variant (`inputs_crescent_nested`) enables CRESCENT's internal AMR refinement (`crescent.amr.max_level = 1`,
`crescent.use_crescent_algorithm = 1` and the child box `crescent.L1_box_lo`/`crescent.L1_box_hi`)
while the ERF coupling stays on CRESCENT's base level.

### Quick smoke test

To confirm the coupled build works end-to-end in a few seconds (3 steps, plotfiles off):

```bash
cd models/ERF/Exec/DevTests/ERF_CRESCENT_Test
mpiexec -np 1 /path/to/build/superbuild/coupler-build/windscape \
    inputs_erf max_step=3 erf.plot_int_1=-1
```

A healthy run prints, each coarse step, a coupling report like

```
[Coupling] Transferred waves from CRESCENT (averaged 5×5):
  eta range: [-1.139601771, 1.13959581] m
  u range:   [-1.245107048, 1.280574333] m/s
  v range:   [0, 0] m/s
  w range:   [-1.292194754, 1.268708673] m/s
```

with finite ranges, and ends with `Coarse STEP 3 ends` followed by a clean
`AMReX ... finalized`. (The full reference case is `max_step = 2000`, `stop_time = 10 s`.)

## Standalone executables

Because ERF is kept coupling-agnostic and CRESCENT retains its own driver, three executables build
from one source tree: **ERF-only**, **CRESCENT-only**, and the **coupled** `windscape`. "Coupling
on or off" is simply a matter of which executable you run. After the superbuild:

```bash
cmake --build build/superbuild/crescent-build              # -> crescent-build/Exec/crescent (+ CTest tests)
cmake --build build/superbuild/erf-build --target erf_exec  # -> erf-build/Exec/erf_exec
```

The ERF-only executable is built with `ERF_ENABLE_TERRAIN_VELOCITY=ON`, like the coupled one; it
runs, for example, the wind-over-waves cases in `models/ERF/Exec/CanonicalTests/WindOverWaves`
(see the README there).

## Components and versions

| Component | Path | Upstream | Pinned version |
|---|---|---|---|
| ERF | `models/ERF` | [github.com/erf-model/ERF](https://github.com/erf-model/ERF) (`development`) | upstream `3b886ecc` (PR #3198) + fork branch `wave-public` @ `08705416` (synced 2026-10-03; `tli-fluid/ERF-mrigark`, copy at `windscape-model/ERF-WAVE` `main`) |
| CRESCENT | `models/CRESCENT` | [github.com/windscape-model/CRESCENT](https://github.com/windscape-model/CRESCENT) | CRESCENT 1.1 release tree (`release/v1.1-prep` @ `b34df45`, from standalone `main` @ `faac11d`; synced 2026-10-05) |
| AMReX | `shared/amrex` | [github.com/AMReX-Codes/amrex](https://github.com/AMReX-Codes/amrex) | `25.12` (+ 2 local patches) |

`models/ERF/UPSTREAM.md` records the upstream base commit of ERF, the layers on top of it,
and the recommended units for contributing WindSCAPE capability back to upstream ERF.
`models/ERF` holds ERF-WAVE, the windscape-model fork of ERF with the moving wave lower boundary.

The two local AMReX patches (a per-instance input namespace and a per-instance physical geometry)
are confined to the mesh layer and let two `AmrCore` solvers run as independent peers in one
process.

## Status and known limitations

- **Single box on the coupled level.** The transfer kernels currently require CRESCENT's coupled
  (base) level to be a single box (`crescent.amr.max_grid_size >= n_cell`). Production multi-rank
  runs use ERF's default horizontal decomposition; a layout-agnostic transfer is planned.
- **Wave nesting under coupling** runs at CRESCENT's base level only; the refined level is internal
  to the wave solver and not yet coupled to ERF (a level-aware exchange API is scaffolded). Since
  the 2026-09-04 CRESCENT sync, nested CRESCENT runs (`crescent.amr.max_level = 1`) require exactly
  one MPI rank, so a coupled nested run is serial too.
- **Two-way feedback** is validated for stability and physical plausibility; quantitative
  energy/momentum-budget validation is future work.

## License

WINDSCAPE is distributed under the BSD 3-Clause License; see [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE). Each vendored component retains its own license: AMReX, ERF and CRESCENT
are BSD-3-Clause, and CRESCENT version 1.0 is released under LLNL release number
LLNL-CODE-2024694. WINDSCAPE was developed at Lawrence Livermore National Laboratory.
Version 1.0 is released under LLNL release number LLNL-CODE-2024697; later versions
inherit this license.

## Acknowledgement

This work was supported by the Laboratory Directed Research and Development (LDRD) program at
Lawrence Livermore National Laboratory (project 25-LW-087). Lawrence Livermore National Laboratory
is operated by Lawrence Livermore National Security, LLC, for the U.S. Department of Energy,
National Nuclear Security Administration under Contract DE-AC52-07NA27344.

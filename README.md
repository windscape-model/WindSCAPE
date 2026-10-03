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

**in process, inside a single executable.** A top-level *Coupler* owns the master time loop and
advances ERF and CRESCENT as equal-weight peers that share one build of AMReX. The coupling is
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
  models/ERF/             vendored ERF (coupling-agnostic; one thin hook added)
  models/CRESCENT/          vendored CRESCENT (full solver, including nesting)
  shared/amrex/           one vendored AMReX, shared by both models (two local patches)
  tests/                  standalone unit tests (e.g. per-instance geometry)
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

Key coupling controls in `inputs_erf`:

```
erf.terrain_type         = MovingFittedMesh   # wave-driven, terrain-following lower boundary
prob.custom_terrain_type = "CRESCENTWaves"
erf.crescent_coupling      = true
erf.crescent_input_file    = "inputs_crescent"
erf.crescent_dt_ratio      = 5                  # CRESCENT sub-steps per ERF coarse step
erf.crescent_grid_ratio    = 5                  # ERF:CRESCENT horizontal grid ratio
```

The wave solver reads its own parameters under the `crescent.` prefix from `inputs_crescent`. A nested
variant (`inputs_crescent_nested`) enables CRESCENT's internal AMR refinement (`crescent.amr.max_level=1`
plus an `crescent.L1_box`) while the ERF coupling stays on CRESCENT's base level.

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
    eta range: [-1.14, 1.14] m   u range: [-1.25, 1.28] m/s   w range: [-1.29, 1.27] m/s
```

with finite `eta`/`u`/`w` ranges, and ends with `Coarse STEP 3 ends` followed by a clean
`AMReX ... finalized`. (The full reference case is `max_step = 2000`, `stop_time = 10 s`.)

## Standalone executables

Because ERF is kept coupling-agnostic and CRESCENT retains its own driver, three executables build
from one source tree: **ERF-only**, **CRESCENT-only**, and the **coupled** `windscape`. "Coupling
on or off" is simply a matter of which executable you run.

## Components and versions

| Component | Path | Upstream | Pinned version |
|---|---|---|---|
| ERF | `models/ERF` | [github.com/erf-model/ERF](https://github.com/erf-model/ERF) (`development`) | `7ed8202` (~PR #3197) |
| CRESCENT | `models/CRESCENT` | CRESCENT HOS wave solver | vendored in full |
| AMReX | `shared/amrex` | [github.com/AMReX-Codes/amrex](https://github.com/AMReX-Codes/amrex) | `25.12` (+ 2 local patches) |

The two local AMReX patches (a per-instance input namespace and a per-instance physical geometry)
are confined to the mesh layer and let two `AmrCore` solvers run as independent peers in one
process.

## Status and known limitations

- **Single box on the coupled level.** The transfer kernels currently require CRESCENT's coupled
  (base) level to be a single box (`crescent.amr.max_grid_size >= n_cell`). Production multi-rank
  runs use ERF's default horizontal decomposition; a layout-agnostic transfer is planned.
- **Wave nesting under coupling** runs at CRESCENT's base level only; the refined level is internal
  to the wave solver and not yet coupled to ERF (a level-aware exchange API is scaffolded).
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
National Nuclear Security Administration under Contract DE-AC52-07NA27344. (The `LLNL-CODE` release
identifier and final auspices wording are to be confirmed with LLNL's software-release process.)

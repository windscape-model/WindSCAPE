# Upstream provenance of the vendored ERF (and how to upstream our work)

WindSCAPE's `models/ERF` is a fork of the public ERF repository. This file is the durable
record of exactly where it forked, what has been layered on top, and how to carry WindSCAPE
capability back to the public repository later.

## Upstream base (the key fact)

| | |
|---|---|
| Upstream repository | <https://github.com/erf-model/ERF> (BSD-3-Clause) |
| Upstream branch | `development` |
| Base commit | `3b886eccd5450aadce06eece455d99949566e13a` |
| Base context | "Refactor Kessler microphysics into BFB-validated helpers (#3198)", upstream state of 2026-05-22 |

Every line of WindSCAPE's ERF is either upstream `3b886ecc` or one of the two layers below, so
a diff of `models/ERF` against upstream `3b886ecc` (minus the vendoring exclusions listed
below) is precisely "our" contribution. Records before 2026-10-03 named `7ed8202` (#3197) as
the fork point; the vendored states of 2026-07-17 to 2026-09-25 already contained `3b886ecc`,
one upstream commit later.

## Layers on top of the upstream base

1. **The wave code of the fork** (since 2026-10-03): fork `tli-fluid/ERF-mrigark`, branch
   `wave-public` (also hosted as `windscape-model/ERF-WAVE`, branch `main`), commit `08705416120ee12f87a1318179ce32d02866178a`, 10 commits on `3b886ecc`:
   `38d03b29` the code of the JCP manuscript (Source/ and the build files of the fork tag
   `jcp-manuscript-2026-10-02` = `eff4101c` on `jcp-writing-pass`); `ffc23fc5` options that
   no published result uses removed (`erf.substep_schedule`, the alias
   `erf.substep_time_order`, `erf.moving_geometry_update = PerStep`,
   `erf.use_acoustic_radiation_top`, a print-only restart block); `c82ef8b7` simplifications
   without a change in results (one wave definition, `Source/Utils/ERF_MovingWave.H`, two
   wave-drag stress functions, one build flag); `5635ad7c` comments only; `a5eb20ae` the
   `MovingSineWave*` surfaces use `prob_lo` and `erf.most.wave_drag` requires
   `MovingSineWaveRamp`; `124fc77f` the wind-over-waves decks in
   `Exec/CanonicalTests/WindOverWaves`; `a9afc0f6` upstream blank lines kept; `fbc4a5f4` the
   orbital velocity written out in each wave-drag stress function; `08705416` a README note on
   what the fork adds; `2fd5703b` documentation
   only (the wind-over-waves README with the four local verification cases, and manual
   entries for the moving-mesh inputs, `erf.beta_d`, `erf.face_sync_int`, the wave drag and
   the terrain-velocity build switch).
   Main content: the moving-mesh (MovingFittedMesh) time integration with the acoustic
   sub-step correction and the GCL fixes, the moving-mesh options that remain
   (`erf.substep_metric_sampling`, `erf.grid_velocity_sampling`, `erf.slow_geometry_time`,
   ...), the orbital moving-wall boundary condition (gated by `ERF_ENABLE_TERRAIN_VELOCITY` /
   `ERF_USE_TERRAIN_VELOCITY`), the wave-drag surface layer, and general corrections that also
   change other ERF cases (listed in the fork's `code_review_wave_2026-10-02.md`, item R4).
   The research material of `jcp-writing-pass` (study folders in `Exec/`, `Docs/ITEM3_*`,
   `mrigark_paper/`, `CLAUDE.md`) is not in this branch; the upstream `paper/` is restored.
   Results are bitwise identical to the tag on CPU and GPU (fork record
   `cleanup_bitwise_2026-10-02/README.md`). The vendored states up to 2026-09-25 came from
   `jcp-writing-pass` (`6717d76`, `2a7c63b`, `ba75feea`).

2. **WindSCAPE coupling hooks** (thin, coupling-agnostic; ERF builds and runs standalone
   with them present). Complete list of touched files:
   - `Source/ERF.H`: `setCouplingHook` + coupling accessor declarations, `m_coupling_hook`
     and wave-surface members (three whole-domain surface slabs: previous, old, new coarse
     step); `setWaveSurfaceVelocity` and the wave-surface x-velocity slabs for the orbital
     moving-wall BC.
   - `Source/ERF.cpp`: one hook invocation in `Evolve()` (after `ComputeDt`, before
     `timeStep`), plus `extractSurfacePressure` / `setWaveElevation` / `fillWaveTerrain` /
     `setWaveSurfaceVelocity` implementations. `fillWaveTerrain` interpolates linearly on
     (old, new) inside the coarse step and on (previous, old) for the lagged surface that
     `erf.grid_velocity_sampling = interval_start` requests.
   - `Source/TimeIntegration/ERF_TI_utils.H`: `fillWaveTerrain` override of every terrain
     surface built from `prob->init_terrain_surface`: the old/src/new (and lagged) surfaces
     of `update_terrain_stage` (MovingFittedMesh per-RK-stage update; `src` at `src_time`,
     which follows `erf.slow_geometry_time`) and the old/new (and previous) surfaces of
     `update_terrain_substep` (acoustic sub-steps; reached only with substepping).
   - `Source/Initialization/ERF_InitCustomTerrain.cpp`: `"CRESCENTWaves"` custom terrain.
   - `Source/BoundaryConditions/ERF_WaveSurfaceVelBC.H` (new): view struct for the
     externally supplied wave-surface x-velocity.
   - `Source/BoundaryConditions/ERF_PhysBCFunct.H` + `Source/ERF_MakeNewArrays.cpp`:
     `ERFPhysBCFunct_u` carries an optional pointer to that view (defaulted, so upstream
     call sites are unaffected), next to the `MovingWave m_wave` of the analytic wall velocity.
   - `Source/BoundaryConditions/ERF_BoundaryConditionsXvel.cpp`: inside the fork's gated
     `ERF_USE_TERRAIN_VELOCITY` branch, prefer the wave-model surface velocity
     (time-interpolated) over the analytic `MovingWave` velocity when the view is active.
   - `Exec/DevTests/ERF_CRESCENT_Test/`: the coupled reference problem (provides
     `amrex_probinit` for the coupled executable), including the `inputs_erf_orbital`
     no-slip moving-wall variant.

## Suggested upstreaming units

When contributing back to `erf-model/ERF`, split the delta into independent PRs:

1. General corrections: the fork keeps them as ten one-commit branches `upstream/*` on
   upstream `development` (fork record `cleanup_bitwise_2026-10-02/README.md`).
2. The moving-wave code (`wave-public`): moving-mesh time integration, wave drag, orbital
   moving-wall BC and the `WindOverWaves` cases, as a feature PR.
3. The generic coupling hook API (`setCouplingHook`, `extractSurfacePressure`,
   `setWaveElevation`/`fillWaveTerrain`, `setWaveSurfaceVelocity` + the
   `WaveSurfaceVelBC` plumbing into the u-velocity BC) — small API PR that would let
   WindSCAPE consume stock ERF; independent of CRESCENT specifics (only the
   `"CRESCENTWaves"` terrain string and the test problem are WindSCAPE-specific).
   Natural companion to unit 2, which it extends from analytic to wave-model data.

For each unit: `git diff 3b886ecc -- <paths>` inside a checkout of upstream with this tree
overlaid, or cherry-pick from the fork branches, whose history is preserved at
`tli-fluid/ERF-mrigark` (the `wave-public` history also at `windscape-model/ERF-WAVE`).

## Vendoring exclusions (not in this tree)

Dropped for repository size, unchanged from upstream:

- `Tests/ERFGoldFiles/` (~428 MB regression baselines) — re-fetch from upstream
  `erf-model/ERF@3b886ecc` (or newer) to run the ERF regression suite.
- `Exec/**/MetGrid/*.nc` real-data inputs (~33 MB) and
  `.Exec_dev/.../ERF_terrain_file.txt` (~11 MB) — re-fetch from upstream.
- The submodules under `Submodules/` (WindSCAPE uses `shared/amrex`, see below).

`mrigark_paper/` and the fork's `CLAUDE.md` are no longer in the source branch.

Force-added although an ignore rule hides them (all four are upstream files):
`Docs/sphinx_doc/Makefile` and `.Exec_dev/DataAssimilation/ReadPltFileAndInterpolate/Makefile`
(ERF's own `.gitignore` rule `Makefile`), and `Build/cmake_with_particles.sh` and
`Build/Perlmutter/cmake_with_cuda_perlmutter.sh~` (the superproject rule `**/build/`, which
also matches `Build/` because git on macOS ignores case, `core.ignorecase = true`). The
vendored states before 2026-10-03 lacked them. After copying a new tree, check
`git status --ignored models/ERF` for further hidden files.

## AMReX note

`wave-public` keeps upstream's `Submodules/AMReX` pin `6f2c4aff90ae3b8a5df8877008de4bc21a71b00b`;
WindSCAPE does not use that pin. All WindSCAPE builds use the single shared `shared/amrex`
(25.12 + the two mesh-layer patches described in the top-level README/NOTICE).

## How this tree was produced (recipe for the next sync)

2026-10-03, from `jcp-writing-pass @ ba75feea + hooks` to `wave-public @ fbc4a5f4 + hooks`.
`wave-public` does not descend from `ba75feea`, so a merge would have brought back removed
options and kept the research files. Instead:

1. In a scratch clone of the fork, branch from `ba75feea`, copy the vendored `models/ERF` over
   it (`rsync --delete`; the vendoring exclusions and `Submodules/` protected) and commit. The
   diff against `ba75feea` is exactly the hook files above and the two marker files.
2. Branch from `wave-public` and apply that delta with `git apply --3way`. Conflicts, all
   resolved on the `wave-public` code:
   - `ERF_TI_utils.H`: `erf.moving_geometry_update` is gone and the `src` surface time is
     `src_time` (`erf.slow_geometry_time`); the override of old/src/new (+ lagged) and of
     old/new (+ previous) is kept; no per-step path.
   - `ERF_PhysBCFunct.H`: the view pointer is added next to `MovingWave m_wave`, which the
     constructor still fills from the inputs under `ERF_USE_TERRAIN_VELOCITY`.
   - `ERF_BoundaryConditionsXvel.cpp`: the wave-model branch wraps the analytic `MovingWave`
     velocity, which stays unchanged.
   `ERF.H`, `ERF.cpp`, `ERF_MakeNewArrays.cpp` and `ERF_InitCustomTerrain.cpp` applied without
   conflict.
3. Copy the result back with the same exclusions, force-add the hidden files, refresh these
   markers, rebuild, and run Tier 0 and the 20-step coupled comparison (`Documentation/CI.md`).

Later on 2026-10-03, the documentation commits `2fd5703b` and `08705416` were applied to
`models/ERF` as the diffs `fbc4a5f4..2fd5703b` and `2fd5703b..08705416` (README.rst); they touch
no hook file and no compiled file.

Next sync from `wave-public`: the vendored state descends from `08705416` again, so the
ordinary three-way merge works: in a scratch clone, branch from `08705416`, overlay
`models/ERF` (records the WindSCAPE delta as one commit), `git merge <new wave-public tip>`,
resolve, copy back, and verify.

# WindSCAPE CI design

Status (2026-10-03): the repository lives at github.com/windscape-model/WindSCAPE, and
the **active CI is GitHub Actions**: `.github/workflows/ci.yml` builds the full superbuild
on ubuntu-latest (ccache-accelerated) and runs the complete Tier 0 smoke suite on every push
to `main` and to any `feature/sync-standalone-*` branch, every pull request, and on manual
dispatch. Every Tier 0 recipe was verified locally on the merged tree before activation
(2026-07-17) and again on the 2026-09-04 standalone sync, which added CRESCENT's CTest
suite and the coupled sub-stepping smoke to Tier 0. The 2026-09-25 sync and the two
2026-10-03 syncs (CRESCENT, then ERF) passed the same local Tier 0 run (8/8).

An **LC GitLab pipeline is also implemented but dormant** (root `.gitlab-ci.yml` plus
`.gitlab/LC/`: pipelines.yml, runners/dane.yml, jobs/dane.yml, `windscape_test.sh`),
following the ERF LC layout from the LLNL `llnl/development` line (parent pipeline
triggering per-machine child pipelines, Jacamar salloc/srun allocation, `lc-templates/
id_tokens` include, `ALLOC_BANK: wwinter`). It activates automatically if this repository
is ever hosted or mirrored on lc.llnl.gov GitLab; GitHub ignores these files. Dane's
Jacamar runners are bound to LC GitLab, so there is no supported way to drive dane jobs
from GitHub; what works without any GitLab is running the same test script manually:

```bash
# on dane, inside an allocation (bank wwinter):
git clone git@github.com:windscape-model/WindSCAPE.git && cd WindSCAPE
salloc -N 1 -t 60 -A wwinter
MODULE_LIST="cmake/3.30.5 ninja/1.11.1 gcc/12.1.1 fftw/3.3.10" \
    bash .gitlab/LC/windscape_test.sh
```

(the script defaults `MPIEXEC_EXECUTABLE=srun` and tolerates running outside a pipeline).

## What CI must protect

WindSCAPE ships three executables from one tree, so CI has three build products and two
behavioral surfaces to guard:

1. The **superbuild** itself: shared AMReX install, `erf_api`, `crescent_static` (+
   standalone `crescent` exe), and the coupled `windscape` link.
2. **Coupled behavior**: the ERF <-> CRESCENT exchange and the moving-terrain feedback.
3. **Component behavior**: each model's own physics, using each model's own test
   infrastructure (both upstreams already have CI we can reuse; see below).

## Tier 0: build + smoke (every push, ~minutes after a warm AMReX cache)

All of these pass locally on the merged tree (first run 2026-07-17; last run 2026-10-03, after
the ERF sync).

| Job | Recipe | Pass criterion |
|---|---|---|
| superbuild | `cmake -S . -B build/superbuild && cmake --build build/superbuild -j` | exit 0; three artifacts exist |
| coupled smoke | in `models/ERF/Exec/DevTests/ERF_CRESCENT_Test`: `mpiexec -np 1 windscape inputs_erf max_step=3 erf.plot_int_1=-1` | exit 0; three `[Coupling] Transferred waves` reports with finite ranges; `Coarse STEP 3 ends` |
| coupled nested smoke | same, plus `erf.crescent_input_file=inputs_crescent_nested max_step=5` | exit 0; five transfer reports; clean finalize (exercises CRESCENT AMR nesting under coupling) |
| CRESCENT standalone smoke | `crescent test_case/inputs_stokes_test max_step=3 crescent.plot_int=-1 crescent.tecplot_int=-1` (mirrors the standalone repo's GitHub CI) | exit 0; `Done timestep 3` |
| orbital moving-wall smoke | `windscape inputs_erf_orbital max_step=20 erf.plot_int_1=-1` (no-slip wall + viscosity; exercises the wave orbital-velocity BC) | exit 0; `Coarse STEP 20 ends` |
| coupled sub-stepping smoke | `windscape inputs_erf max_step=5 erf.plot_int_1=-1 erf.substepping_type=Implicit erf.fixed_mri_dt_ratio=2 erf.use_lagged_delta_rt=false` (exercises the wave-surface override of the acoustic sub-step terrain update, which the default deck never reaches) | exit 0; `Coarse STEP 5 ends` |
| CRESCENT CTest suite | `cmake --build build/superbuild/crescent-build && ctest --test-dir build/superbuild/crescent-build` (the standalone repo's unit tests and input-contract rejection decks, 22 cases, built against the shared AMReX) | all tests pass |
| per-instance geometry unit test | build + run `tests/per_instance_geom` against the shared AMReX | its own asserts (two AmrCore instances with independent domains/periodicity across AMR levels) |

Cache note: the shared AMReX build dominates cold-start time. Cache
`build/superbuild/install/amrex` keyed on `shared/amrex` content (or use ccache); ERF is the
next-largest unit.

## Tier 1: coupled regression gate (ACTIVE, in the GitHub Actions job)

Runs in the **same ubuntu job as the Tier 0 smokes**, right after them: a 20-step coupled
run with a plotfile at step 20, compared with the superbuild's own
`install/amrex/bin/amrex_fcompare` against the committed golden:

```
mpiexec -np 1 windscape inputs_erf max_step=20 erf.plot_int_1=20
amrex_fcompare -r 1e-12 --abs_tol 1e-13 tests/golden/linux-x86_64/plt00020 plt00020
```

- **Where the golden comes from (and why it lives where it does).** Bitwise-level goldens
  are platform + compiler specific, so the golden is generated **on the CI runner itself**
  (linux x86_64, distro gcc), not on a developer machine: trigger the workflow manually
  with the `regenerate_golden` input (`gh workflow run CI -f regenerate_golden=true`).
  The dispatch run commits the fresh golden to the branch itself (workflow token,
  `permissions: contents: write`) as `tests/golden/linux-x86_64/plt00020` (~3 MB); that
  push then triggers a normal CI run that exercises the new golden. The gate step
  activates automatically whenever the golden directory exists and is skipped on
  regeneration dispatches.
- Refresh policy: the golden changes only when physics intentionally changes (for example a
  standalone sync like 2026-07-17); refresh it in the same change and say so in its message.
  The 2026-09-04 sync (CRESCENT f65ac7f -> f228bf6, ERF fork 6717d76 -> 2a7c63b) left the
  reference case unchanged to round-off: locally, plt00000 is bitwise identical to the
  pre-sync binary and step 20 agrees within 1.4e-12 relative (absolute <= 3e-10), and on
  the runner the gate passed against the 2026-07-17 golden without regeneration (run
  33938600232, first push of branch feature/sync-standalone-20260904; ~22 min with the
  July ccache). The 2026-10-03 sync (CRESCENT f39cacf -> f15a84e, the 2026-10-02 code
  review; ERF fork unchanged) left the reference case bitwise unchanged: locally, plt00000
  and plt00020 equal the pre-sync binary's in every variable, and the logs of the four
  coupled smokes differ only in wall time and three message texts. The 2026-10-03 ERF sync
  (fork `jcp-writing-pass` ba75feea -> `wave-public` fbc4a5f4, the cleaned manuscript code with
  the WindSCAPE hooks re-applied) also left the reference case bitwise unchanged on macOS, and
  with it the orbital, sub-stepping and `IntervalStart` variants; the smoke logs are identical.
  On the runner the Tier 1 gate passed against the unchanged golden after both 2026-10-03
  syncs (runs 37129426844 and 37143386976). The documentation updates of the same day
  (`wave-public` 2fd5703b and 08705416) change no compiled file.
- Tolerances: same-platform reruns of identical binaries reproduce bitwise, so
  `-r 1e-12 --abs_tol 1e-13` is effectively exact while absorbing pure last-digit noise
  (for example a runner-image compiler bump). Real changes sit far above it: the 2026-07-17
  sync moved step-20 velocities by up to 4e-3 relative. On a gate failure the candidate
  plotfile is uploaded as the `regression-candidate-plt00020` artifact for offline fcompare.
- macOS developers should not compare against this golden bitwise; run the same case
  locally against a locally generated reference, or use a loose tolerance (~1e-10).
  The dormant LC pipeline retains the separate gold-repo clone pattern for LC machines.

## Tier 2: reusing each model's own CI inside WindSCAPE

### ERF (answers "can we use ERF's CI even though ERF is coupled with CRESCENT?")

Yes. The vendored ERF is deliberately **coupling-agnostic**: the hooks are inert unless a
driver registers them, and `models/ERF` still configures/builds/tests with its own CMake.
Two adaptations are required because WindSCAPE vendors sources and shares one AMReX:

1. Build ERF's test suite against the shared AMReX:
   `cmake -S models/ERF -B <bld> -DERF_USE_INTERNAL_AMREX=OFF -DAMReX_DIR=<install>/lib/cmake/AMReX -DERF_ENABLE_TESTS=ON ... && ctest --test-dir <bld>`
   (drop `GIT_SUBMODULE_STRATEGY: recursive` from any adapted CI YAML; WindSCAPE has no
   submodules, and `amrex_fcompare` comes from the superbuild install, not
   `Submodules/AMReX`).
2. Gold files: WindSCAPE dropped `Tests/ERFGoldFiles` (~428 MB) from the tree. Upstream's LC
   CI already solves this by cloning a gold repo and setting
   `ERF_TEST_GOLD_FILES_DIRECTORY`; do the same. Important nuance: our ERF carries the
   fork's moving-mesh and wave-drag changes and its general corrections
   (`models/ERF/UPSTREAM.md`), so upstream's gold answers may legitimately differ for some of
   its regression tests. Recommended: generate a
   WindSCAPE-owned gold set once from the merged tree (same machine/compiler policy as
   upstream's per-machine gold dirs) and treat upstream gold as an advisory comparison.

The complete LC GitLab machinery is already vendored at `models/ERF/.gitlab/LC/`
(`.gitlab-ci.yml`, `gitlab_test.sh`, runner/job YAMLs for dane, tioga, tuolumne, matrix). It
is upstream's, unmodified by the fork, and is the natural template for WindSCAPE's own LC
pipeline; the project-specific pieces (runner tags, `ALLOC_BANK`, `ALLOC_QUEUE`, module
lists, gold-repo URL/token grants) must come from our LLNL GitLab project settings.

### CRESCENT

The standalone repo's CI (vendored at `models/CRESCENT/.github/workflows/ci.yml`) is a
build + 3-step smoke; Tier 0 already reproduces it against the shared AMReX, and since
2026-09-04 Tier 0 also runs the standalone repo's CTest suite (`models/CRESCENT/tests/`,
enabled by `include(CTest)` in CRESCENT's CMake; the superbuild builds the test executables
with `cmake --build build/superbuild/crescent-build`). CRESCENT's bit-for-bit regression gate
and its physics validation decks live in the standalone development repository, not in
WindSCAPE, which vendors the CRESCENT release tree (since 2026-10-03); the gate's baselines
(~177 MB, deterministically regenerable) are too heavy for per-push CI.

## CI file map

| File | Status | Role |
|---|---|---|
| `.github/workflows/ci.yml` | **active** | GitHub Actions: superbuild + full Tier 0 smoke suite (coupled, nested, orbital, sub-stepping, standalone CRESCENT, CRESCENT CTest suite, per-instance geometry) + the Tier 1 golden gate, ccache-cached |
| `.gitlab-ci.yml` (root) | dormant | LC parent pipeline: `ALLOC_*` variables (bank `wwinter`), workflow rules, per-machine trigger template |
| `.gitlab/LC/pipelines.yml` | dormant | One child-pipeline entry per machine (dane; copy ERF's runner/job YAMLs to add tioga/tuolumne) |
| `.gitlab/LC/runners/dane.yml` | dormant | Adopted from ERF `llnl/development`: machine-up check, `salloc` named allocation, `srun` of the test script, release stage |
| `.gitlab/LC/jobs/dane.yml` | dormant | Toolchain job (`cmake/3.30.5 ninja/1.11.1 gcc/12.1.1 fftw/3.3.10`) |
| `.gitlab/LC/windscape_test.sh` | usable manually on LC | Builds the superbuild, runs all Tier 0 smokes and, when `WINDSCAPE_GOLD_REPO` is set, the Tier 1 fcompare gate |

Tier 2 (ERF's own CTest suite with an external gold clone) remains a design item; add it as
a scheduled job once a gold-files location exists.

## Open items

- If the LC pipeline is ever activated (repo mirrored to lc.llnl.gov): verify the dane FFTW
  module name in `.gitlab/LC/jobs/dane.yml` on the first run (the queue/bank are set:
  `pci` / `wwinter`).
- Decide the cadence for the CRESCENT regeneration gate and the ERF full CTest job
  (nightly vs weekly), and whether they run on Actions (portable parts) or manually on LC.

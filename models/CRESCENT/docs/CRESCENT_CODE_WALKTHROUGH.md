# CRESCENT code walkthrough

A guided tour of the source tree for developers: how the solver is laid
out, how a time step flows through the code, where each piece of the
algorithm lives, and how to extend it safely. Companion documents:

| Document | What it covers |
|---|---|
| `CRESCENT_ALGORITHM_OVERVIEW.md` | The algorithm itself, method-level |
| `R2_CORRECTION_FIELD_SPEC.md` | The CCF coupling, equation-level spec |
| this file | The code: files, call paths, design rules |

Paper-to-code vocabulary: the paper's *blending function* sigma is
`mask`/`sigma` in the code, the *attenuation function* chi is the
complementary `taper`/`chi`, the *accurate region* Omega_A is `roi` /
"accurate" in identifiers, and the *auxiliary copies* are the periodic
`scratch` fields. The code predates the paper's register audit and
keeps its own names; do not rename code to chase the paper.

## 1. The thirty-second tour

One executable, `build/Exec/crescent`, driven by a single ParmParse input
file (`crescent.*` keys). Two levels: a doubly periodic parent L0 solved
by a classical HOS pseudo-spectral method, and an optional one-way nested
child L1 in a rectangular window. Each step: L0 advances four RK4 stages;
L1 then advances its own four stages, coupling to recorded parent-stage
data; a post-step finalization anchors the child to the parent; outputs
and per-step diagnostic CSVs are written on their cadences.

```
main.cpp
  └─ CRESCENT (Source/CRESCENT.H, an amrex::AmrCore)
      ├─ Init()                     level creation, wave init, workspaces
      └─ Run()                      solver/CRESCENT_RunLoop.cpp
           └─ per step: advanceAllLevelsForStep  (solver/CRESCENT_TimestepSupport.cpp)
                ├─ L0: core::free_surface_evolution   4 x runRKStage
                ├─ L1: nested::free_surface_evolution::advanceL1Timestep
                │        4 x (RK4 substep + total_field::runStage
                │             or correction_field::runL1CorrectionFieldStage)
                │        RK4 combination + state relaxation or correction taper
                ├─ outputs: plot / tecplot / movie-tecplot / checkpoint
                └─ diagnostics: energy, oscillation, edge-probe CSVs
```

## 2. Directory map

| Directory | Role |
|---|---|
| `Source/main.cpp`, `Source/CRESCENT.H` | Entry point and the solver facade: an `amrex::AmrCore` subclass owning levels, state, and the AMR callbacks (`MakeNewLevelFromScratch`, `MakeNewLevelFromCoarse`, ...) |
| `Source/core/` | The periodic L0 HOS solver: RK4 dynamics, spectral operators, order-M surface vertical-velocity series (legacy `DNO` source names) |
| `Source/nested/` | Everything L1: window geometry, masks, both couplings, per-stage orchestration, nested diagnostics |
| `Source/solver/` | The driver: run loop, setup, inputs, checkpoint/restart, outputs, step diagnostics |
| `Source/shared/` | Cross-cutting infrastructure: solver state and the two contexts, solution state, spectral workspace, weight tables and kernels, interpolation, I/O |
| `Source/init/` | Wave initialization: Stokes, JONSWAP, crescent seed |
| `Source/external/` | A surface-exchange interface for an external coupler (not used by the solver itself) |
| `Source/archive/` | Dead prototypes, excluded from the build; ignore |
| `tests/` | CTest unit tests and invalid-input decks (built when `BUILD_TESTING` is on) |
| `tools/movie/` | The Tecplot-to-movie pipeline |

## 3. Execution path, file by file

**Startup.** `main.cpp` checks the input file exists, then calls
`amrex::Initialize` with `applyAmrDefaultInputs`, which pre-loads AMReX
ParmParse defaults (grid layout, `max_grid_size`, `blocking_factor`) so
decks stay small. It constructs one `CRESCENT` object and calls `Init()`
then `Run()`.

**Setup.** The `CRESCENT` constructor calls `constructSolverState`
(`solver/CRESCENT_Setup.cpp`), which sizes the per-level containers and
reads every input through `readParameters`
(`solver/CRESCENT_SetupParameters.cpp`: run control, output, checkpoint
and wave keys) and `nested/CRESCENT_NestedAlgorithmInputs.cpp` (the
nested keys, into `AlgorithmConfig` in `nested/CRESCENT_NestedConfig.H`).
Both validate in place; a deck that sets an invalid or retired value
stops at t = 0. After AMReX builds the hierarchy, `Init()` calls
`initializeSolverState` (fresh run) or the checkpoint reader, which calls
`rebuildRestartDerivedState` (restart); both are in
`solver/CRESCENT_Setup.cpp`. They create the per-level FFT plans
(`amrex::FFT::R2C`, FFTW on hosts and cuFFT under CUDA) and workspaces,
initialize the parent wave field (`init/`), and set up the child
(`nested/CRESCENT_NestedSetup.cpp`: coupling zones, storage, the child
state interpolated from the parent).
`validateAlgorithmSelection` enforces the startup invariants, for example
a single MPI rank for nested runs and, in mode 2, equal state and
derivative interpolation orders.

**The run loop** (`solver/CRESCENT_RunLoop.cpp`). Owns an
`ExecutionState` (time, last plot/tecplot/checkpoint step), writes
initial outputs, then loops: advance one step, honor each output
cadence, emit step diagnostics, checkpoint. Restart re-creates the
output directory and resumes cadences from the checkpointed step
(`solver/CRESCENT_Checkpoint.cpp`). Checkpoint format v10 serializes all active
wave/spectral and nested evolution settings (the slots of retired options
hold their supported value); derived fields are rebuilt by
the coupling-specific path, and checkpointed settings override conflicting
restart deck values.

**One L0 step** (`core/CRESCENT_CoreFreeSurfaceEvolution.cpp`).
`runRKStage` is called four times. Each stage: `dynamics::advanceSurfaceRK4`
(`core/CRESCENT_Dynamics.cpp`) applies the stage update from
`(eta_saved, phi_saved, rhs_eta, rhs_phi)`; `spectral::applyLowPassFilter`
(`core/CRESCENT_Spectral.cpp`) applies the parent low-pass filter;
`dno::computeSolutionStateDerivativeBundle`
(`shared/CRESCENT_DNOExecution.cpp`) computes the surface derivatives and
the order-M surface vertical velocity (`physics::computeVerticalVelocity`
in `core/CRESCENT_Physics.cpp`: eta powers, the potential expansion and
the vertical-velocity series); `dynamics::computeSurfaceRHS` assembles
and dealiases the tendency. Per-mode weights of the filters and the
|k|^n factors of the series come from cached tables
(`shared/CRESCENT_SpectralWeightTable.H`). In CCF runs a hook captures the
parent's per-stage fields for the child
(`nested/CRESCENT_NestedCorrectionField.cpp`, `captureCoreStage`); the
capture does not alter L0. After the step, the parent derivative bundle
of the new state is computed once and kept as the interpolation source of
the child (`nested/CRESCENT_NestedCoreSnapshots.cpp`).

**One L1 step** (`nested/CRESCENT_NestedTimestep.cpp`).
`advanceL1Timestep` runs the four stages with the classical weights
`c = {0, 1/2, 1/2, 1}`. Each stage applies the RK4 update
(`advanceSurfaceRK4`, the same kernel as L0) and then the stage of the
selected coupling (`crescent.l1_coupling_mode`):

```cpp
for (int irk = 0; irk < 4; ++irk) {
    dynamics::advanceSurfaceRK4(..., irk);
    if (correction_field) {
        correction_field::runL1CorrectionFieldStage(
            level_context, irk, time + c[irk] * dt, c[irk]);
    } else {
        total_field::runStage(level_context, irk);
    }
}
dynamics::updateSurface(...);
// then correction_field::applyPostRK4CorrectionFieldTaper (CCF)
// or total_field::applyStateRelaxation (PTF)
```

- **Mode 1, PTF (production default)**, `nested/CRESCENT_NestedTotalField.cpp`.
  In call order: copy the live state to the auxiliary periodic copies;
  regularized C3-target interface repair on those copies
  (`NestedC3InterfaceRepair.cpp`, one regularized least-squares line fit
  per grid line, x pass then y pass; the small dense solver is
  `NestedC3LinearSolvers.H`, the basis bookkeeping `NestedC3Periodicity.H`,
  the eta/phi pair dispatch `NestedC3Dispatch.H`); the mean-preserving
  combined Hou-Li / two-thirds / low-pass filter (`NestedFilter.cpp`);
  native spectral derivatives (`NestedSurfaceDerivatives.cpp`) and the
  full-order surface vertical velocity of the periodic pair
  (`dno::computeNestedL1DNO`); then the tendency blend (`NestedRK4.cpp`),
  combining the child tendency with the stage-exact interpolated parent
  tendency through sigma. Stage 0 reuses the preparation of the previous
  step end, which acted on the same state.
- **Mode 2, CCF**,
  `nested/CRESCENT_NestedCorrectionField.cpp`. No periodic copies and no
  C3 repair exist on this path; every transformed field is compact or a
  Fourier series by construction. The file contains the whole coupling:
  reconstruction of the pre-filter parent stage baseline, the attenuated
  correction pair, child FFT derivatives and the correction
  vertical-velocity map with its order-2 background-correction cross
  terms (spectral utilities such as `applySpectralRieszMultiplier`, the
  |k| half-Laplacian symbol with optional 2/3 cutoff, are local to the
  file), the child-assembled background tendency (the resolution
  commutator that seeds the correction), the optional Dommermuth
  adjustment ramp, the projection onto the added band
  (`applyAddedBandProjector`), and the assembled stage RHS with
  parent-parity dealiasing. Background-containing metric products are
  re-tapered before their child projections, and the final projection of
  the child vertical velocity is bracketed by tapers. Applied pressure is
  kept outside spectral projection, as on the parent and PTF paths.

**Post-step finalization.** CCF re-anchors the state to the
decomposition, interpolated parent plus attenuated correction (the
projection of the paper's Eq. eq:ccf-attenuation); PTF applies the
repaired-and-filtered state relaxation toward the interpolated parent.
The relaxation is the necessary energy drain of the PTF path; removing
it degrades grid-scale noise by an order of magnitude (an A/B result of
the development repository). The end-of-step output bundle of the child is then prepared by
`free_surface_evolution::recomputeL1OutputDerivatives`, the one dispatch
point for both couplings.

## 4. Solver state and contexts

`CRESCENT` owns one `SolverState` (`shared/CRESCENT_SolverState.H`):
`CoreSolverState` (run control, wave parameters, per-level fields and
spectral workspaces) and `NestedSolverState` (the nested configuration and
the child coupling state). Code reaches it through two small contexts in
`shared/CRESCENT_SolverContext.H`:

| Context | Gives |
|---|---|
| `SolverContext` | the whole hierarchy: state, mesh, finest level; `levelContext(lev)` |
| `LevelContext` | one level: `level()`, `core()`, `nested()`, `solutionState()`, `runtime()`, `spectralWorkspace()`, `nestedConfig()`, `nestedRuntime()`, the level and parent meshes and geometries, `waveOperator()` and `coreWaveOperator()` |
| `LevelWaveOperatorView` | the per-level spectral operator bundle (FFT plan, workspace, wavenumbers, filter settings) used inside a stage |

Physics kernels take MultiFabs and plain values, not contexts. The
2026-10 code review replaced seven context classes, six access structs
and their factories by these two contexts; none of the former classes
restricted access, because each exposed the whole core and nested state.

## 5. Data model

- `shared/CRESCENT_SolutionState.H`: per-level `eta`, `phi`, the
  stage-saved pair `eta_saved`/`phi_saved`, and the stage tendencies
  `rhs_eta`/`rhs_phi`, all `amrex::MultiFab`s. The RK4 kernel is written
  against exactly this quintet.
- `shared/CRESCENT_WaveState.H` and `shared/CRESCENT_SpectralWorkspace.H`:
  derived wave quantities and the FFT plans/buffers per level.
- `nested/CRESCENT_NestedTypes.H` and `NestedState.H`: everything the
  child needs beyond the solution: the blending function sigma and its
  complement chi (`NestedMasks.cpp` builds them; C3-smooth septic
  profile, zones Omega_A / Omega_M / Omega_E from
  `roi_accurate_fraction` / `roi_spectral_fraction`), the window
  geometry (`NestedDomains.cpp`), interpolated parent bundles and
  recorded parent-stage captures, and `coupling_derivatives`, the
  per-stage derivative bundles exchanged between the stage orchestrator
  and the couplings (`NestedDerivativeBundles.H`). CCF retains 11 coarse
  fields for each of four RK stages and the final output state; 10 of those
  quantities are interpolated to the child as the active background bundle.
- Interpolation between levels is `shared/CRESCENT_Interpolation.cpp`
  (the solver's own prolongation; the same operator the analysis
  scripts and the paper's I denote).

The periodicity constraint of the paper is a code rule here: fields in
`SolutionState` on L1 are windows of a larger sea and are never
transformed; only the auxiliary periodic copies (PTF) or the attenuated
correction, re-tapered background products, and outputs of child inverse
transforms (CCF) enter `SpectralWorkspace` plans. Interpolated background
fields themselves remain in physical space. Runtime monitoring of the rule
lives in the diagnostics below.

## 6. Inputs and toggles

All keys use the `crescent.` prefix; nested keys are parsed in
`NestedAlgorithmInputs.cpp` into `AlgorithmConfig`, the other keys in
`solver/CRESCENT_SetupParameters.cpp`. The load-bearing ones:

| Key | Consumed by | Meaning |
|---|---|---|
| `use_crescent_algorithm` | driver | enable the nested L1 path |
| `L1_box_lo/hi` | `NestedDomains` | child window in L0 cell-edge coordinates |
| `roi_accurate_fraction`, `roi_spectral_fraction` | `NestedMasks` | zone geometry (Omega_A margin, spectral margin) |
| `l1_coupling_mode` | `NestedTimestep` | 1 = PTF (default), 2 = CCF |
| `state_interpolation_order`, `derivative_interpolation_order` | `Interpolation` | 1 = bilinear, 2 = quadratic (default) parent-to-child prolongation |
| `perturbation_order` | `core/Physics`, vertical-velocity series | HOS order M (crescent decks use 4) |
| `c3_boundary_penalty`, `c3_coefficient_regularization` | C3 files | PTF periodization weights (the buffer widths follow from the zone geometry); never consulted in mode 2 |
| `beta1_filter`, `beta2_filter`, `beta1_filter_l1` | `Spectral`, `NestedFilter` | low-pass onset/cutoff in multiples of k_peak |
| `l1_rolloff_alpha`, `l1_rolloff_power` | `NestedFilter` | child Hou-Li roll-off exp(-alpha (k/k_N)^p) |
| `adjustment_periods` | `Dynamics`, both couplings | Dommermuth nonlinear ramp in peak periods |
| `peak_wavenumber_override` | `WaveState` | explicit filter peak (-1 = the initializer's peak) |
| `write_omega_a_edge_diagnostics` | solver diagnostics | per-step edge-probe CSV |

Retired keys. The 2026-10 code review removed the derivative blend
(mode 0) and the legacy and exploratory branches of fifty switches: the
alternative masks and C3 variants, the compatibility branches of the
2026-07 to 2026-08 fixes, the legacy wave initializers, the post-timestep
low-pass filter, and the unimplemented NetCDF/HDF5 outputs; on 2026-10-03
`ccf_band_restrict` followed (the projection onto the added band is part of
CCF). A deck may still state a retired key with its supported value; any
other value stops
the run with a message that names commit `f39cacf` of the development
repository, where the other branches remain. `shared/CRESCENT_RetiredOptions.H` implements the rule;
checkpoint format 10 keeps the slots of retired keys with their supported
value. `blending_steepness` never had an effect; it is accepted with a
warning because older decks set it.

Development rule (enforced by review, not by code): physics changes land
behind a new `crescent.*` toggle defaulting to the old behavior, get
A/B-tested on identical decks, and only then flip the default in a
separate commit. Once the evidence is in, retire the old branch rather
than keep it as a permanent switch.

## 7. Diagnostics and outputs

| File (run directory) | Producer | Content |
|---|---|---|
| `crescent_energy.csv` | `nested/CRESCENT_NestedEnergy.cpp` via step diagnostics | per-step canonical Zakharov energies; `G psi` is reconstructed pointwise from the derivative/vertical-velocity bundle, and the fidelity metric is `energy_L1_omegaA` vs `energy_L0_under_omegaA` (whole-child energy is blend-zone dominated and not meaningful) |
| `nested_oscillation_diagnostic.csv` | `solver/CRESCENT_RuntimeDiagnostics.cpp` | second differences of eta in and around Omega_A (sharpest artifact detector) plus wrap-junction periodicity metrics of every FFT input copy |
| `omega_a_edge_diagnostic.csv` | `solver/CRESCENT_RuntimeDiagnostics.cpp` (called from `StepDiagnostics.cpp`) | probe lines at the Omega_A x-edges; part of the regression gate, its columns are frozen (add new files, never new columns) |
| `<prefix>_tecplot*.dat`, `<prefix>_movie_tecplot*.dat` | `shared/CRESCENT_IO.cpp` | full-field and movie (eta-only, L0 full + L1 restricted to Omega_A) Tecplot outputs; movie files feed `tools/movie/` |
| plotfiles, checkpoints | `RuntimeOutputs.cpp`, `Checkpoint.cpp` | AMReX-native outputs and restart state; CCF plot derivatives are total fields, while retained RK tendencies are labeled `rhs_*_last_step_k1` |

Output bundles are populated before step-0 files. CCF parent-interpolation
columns therefore contain the actual interpolated parent derivatives, and
the blended child `u/v` diagnostics are formed pointwise without sending the
non-periodic blended fields through an FFT.

The developer trace logging and the research-only diagnostics were removed
in 2026-09 (logging cleanup) and 2026-10 (seven unread CSVs, among them
the Omega_S corner and interface files). Use `BL_PROFILE` regions and an
AMReX TinyProfiler build for timing.

## 8. Initialization

`init/CRESCENT_WaveInit.cpp` implements the Stokes carrier and the
crescent case (fifth-order carrier plus oblique free-wave seed, ported
from the reference `vort.f90::crescent_init` and validated to ten
digits); `init/CRESCENT_JonswapInit.cpp` the seeded JONSWAP spectra
(long-crested and directional). There is one construction: a normalized
directional spreading density evaluated at signed angles, full-circle
directional support, independent phases where both propagation senses
occupy an R2C pair, and the required Hermitian constraints on the
self-conjugate R2C planes (`init/CRESCENT_JonswapSpectrum.H`);
`jonswap_amplitude_scale`, not `wave_steepness`, is its optional amplitude
multiplier. The Stokes (linear Airy) initializer keeps the signed
propagation quadrant of the requested angle. L1 is never initialized
independently: it is filled by interpolation from L0
(`nested/CRESCENT_NestedSetup.cpp`), which is
why an inert coupling reproduces the interpolated parent exactly, the
signature that exposed the 2026-07-10 CCF bug.

## 9. GPU design rules

One codebase, two backends; `cmake -DCRESCENT_ENABLE_CUDA=ON` compiles
the same sources for device. The port's working rules:

- Never capture host pointers in device lambdas. `geom.CellSize()`
  returns host memory; capture `GpuArray` values (or a `GeometryData`)
  instead. This was the recurring crash class of the port.
- All-device physics: per-cell work goes through `amrex::ParallelFor`;
  transforms through `amrex::FFT::R2C` (cuFFT under CUDA).
- I/O staging: Tecplot gathers stage through `The_Pinned_Arena`
  (`shared/CRESCENT_IO.cpp`); device-to-host copies from pageable
  buffers were the tecplot segfault source.
- Dual-path where CPU bit-reproducibility matters: the C3 repair keeps
  the original scalar host path and adds a batched per-line device
  kernel under `AMREX_USE_GPU` (`NestedC3InterfaceRepair.cpp`,
  `c3SeamRepairPassDevice`; shared solver code is annotated
  `AMREX_GPU_HOST_DEVICE` in `NestedC3LinearSolvers.H`).
- Cached weight tables (`SpectralWeightTable`) live on the FFT's spectral
  layout and are read inside the transform callbacks, on host and device.
- Acceptance: the CPU build must stay bit-for-bit (the development
  repository's regression gate); GPU vs CPU agree to machine precision (~1e-15) on
  full nested runs, never bit-identical (FFT algorithms, parallel
  reductions, and FMA contraction reorder rounding).

## 10. Verification workflow

- **Unit and input tests** (`tests/`, built when `BUILD_TESTING` is on, the
  CMake default): the linear dynamics, energy, C3 least squares, checkpoint
  round trip, interpolation, JONSWAP and Stokes constructions, and the
  rejection of invalid or retired inputs. Run them with
  `ctest --test-dir build --output-on-failure`.
- **Sample decks** (`test_case/`): the three shipped inputs exercise the
  single-domain Stokes and JONSWAP paths and the nested initialization.
- The bit-for-bit regression gate, the A/B suite, the physics validation
  suite and the analysis scripts behind the published results are
  maintained in the development repository and are not part of this
  release.

## 11. Extending the solver

**A new physics variant.** Add a key in `NestedAlgorithmInputs.cpp`
(default = old behavior) and a field in `NestedConfig`; branch at the
narrowest point that distinguishes the variants (stage orchestrator or
inside the coupling); A/B on one validation deck; run the gate (must be
bit-for-bit at defaults); flip the default in its own commit and record
the A/B evidence.

**Respect the periodicity constraint.** Before adding any FFT, filter,
or dealiasing call on the child, answer: is this field one of the
auxiliary periodic copies, an attenuated compact field, a Fourier
series assembled from parent modes, or a product of those? If it is the
live state, an interpolated parent field, or anything sigma-blended,
it may only be used pointwise in physical space. Vertical-velocity-map outputs fill the
domain even for compact inputs and must be re-attenuated by chi before
multiplying non-periodic background data. The wrap-junction columns of
`nested_oscillation_diagnostic.csv` will catch violations at run time;
believe them.

**A new diagnostic.** Write a new CSV; never add columns to
`omega_a_edge_diagnostic.csv` (gate-frozen). Follow the pattern of
`solver/CRESCENT_RuntimeDiagnostics.cpp` (compute, reduce, the I/O rank writes).

**A new wave initializer.** Add to `init/`, register in the
`initializeCoreWaveField` dispatch (`init/CRESCENT_WaveInit.cpp`), and
remember L1 always
initializes by interpolation; a child-only initializer is a design
error in a one-way nested solver.

**A new test.** Add a unit test or an invalid-input deck under `tests/`
and register it in `tests/CMakeLists.txt`.

## 12. Suggested reading order

1. `main.cpp`, `CRESCENT.H`, `solver/CRESCENT_RunLoop.cpp`,
   `solver/CRESCENT_TimestepSupport.cpp` — the skeleton.
2. `core/CRESCENT_CoreFreeSurfaceEvolution.cpp` then `core/CRESCENT_Dynamics.cpp`
   and `core/CRESCENT_Spectral.cpp` — a periodic HOS step end to end.
3. `core/CRESCENT_Physics.cpp` — the order-M surface vertical-velocity
   series, the hardest part of the core.
4. `nested/CRESCENT_NestedMasks.cpp` and `NestedDomains.cpp` — the geometry.
5. `nested/CRESCENT_NestedTimestep.cpp` — the child step and the mode
   dispatch.
6. `nested/CRESCENT_NestedTotalField.cpp` with the C3, filter and RK4
   files it calls — the default coupling.
7. `nested/CRESCENT_NestedCorrectionField.cpp` — the self-contained
   CCF coupling.
8. `nested/CRESCENT_NestedEnergy.cpp` — how fidelity is measured.

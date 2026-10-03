# ERF-CRESCENT Coupling Test

A coupled ERF + CRESCENT air-sea test case. Each step, ERF's surface pressure
forces the waves, and the evolving CRESCENT wave surface moves ERF's
terrain-following lower boundary.

To build and run, see "Running a coupled case" in the top-level WINDSCAPE
`README.md`. In short, build the superproject and run from this directory:

```bash
mpiexec -np 1 /path/to/build/superbuild/coupler-build/windscape inputs_erf
```

Files in this directory:

- `inputs_erf`: ERF parameters; enables the coupling and selects the
  `CRESCENTWaves` wave-driven moving terrain.
- `inputs_crescent`, `inputs_crescent_nested`: CRESCENT wave parameters (the
  nested variant turns on CRESCENT's internal AMR).
- `ERF_Prob.H`, `ERF_Prob.cpp`: problem setup.

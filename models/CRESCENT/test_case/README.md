# Test Cases

This directory is the curated in-repo sample set for CRESCENT runs.

## Compile (from the repository root)

```bash
git submodule update --init --recursive
cmake -S . -B build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build --target crescent -j4
```

The executable is `build/Exec/crescent`. See the root README for external
AMReX options.

## Sample cases

- `inputs_stokes_test`: deterministic single-domain Stokes case for quick
  numerical sanity checks.
- `inputs_CRESCENT_L0_only`: single-domain JONSWAP case that exercises the
  core solver path without nested coupling.
- `inputs_nested_jonswap_init`: two-level nested JONSWAP initialization
  example for CRESCENT setup/output checks.

## Run and find the output

The solver writes all output relative to the current working directory
(plotfiles, Tecplot files, and diagnostic CSVs such as
`crescent_energy.csv`). Always run from inside a dedicated output
directory:

```bash
mkdir -p runs/stokes && cd runs/stokes
../../build/Exec/crescent ../../test_case/inputs_stokes_test > run.log 2>&1
```

Afterwards the directory contains `run.log` and, as the deck selects, AMReX
plotfiles (`<plot_file>NNNNN/`), Tecplot ASCII files (`*_tecplot*.dat`), and,
for nested runs, the energy and diagnostic CSV files. Raw run data is never committed to the
repository: everything under `runs/` is gitignored, and all
results are regenerated from these input decks (runs are deterministic
for fixed seeds and rank counts).

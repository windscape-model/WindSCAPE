# Test Cases

This directory is the curated in-repo sample set for CRESCENT runs.

It intentionally keeps only a few representative inputs:

- `inputs_stokes_test`: deterministic single-domain Stokes case for quick
  numerical sanity checks.
- `inputs_CRESCENT_L0_only`: single-domain JONSWAP case that exercises the core
  solver path without nested coupling.
- `inputs_nested_jonswap_init`: two-level nested JONSWAP initialization example
  for CRESCENT setup/output checks.

Typical run:

```bash
cd build/Exec
./crescent ../../test_case/inputs_stokes_test
```

Use the other input files the same way.

# CRESCENT Documentation

*CRESCENT: A block-structured, adaptive-mesh framework for ocean surface wave modeling*

This directory contains public CRESCENT documentation for users, reviewers, and
downstream developers.

## Documentation

- [CRESCENT_ALGORITHM_OVERVIEW.md](CRESCENT_ALGORITHM_OVERVIEW.md): concise overview of the solver, nested algorithm, checkpoint/restart path, and external exchange interface.
- [CRESCENT_CODE_WALKTHROUGH.md](CRESCENT_CODE_WALKTHROUGH.md): source layout, execution paths, controls, and verification workflow.
- [R2_CORRECTION_FIELD_SPEC.md](R2_CORRECTION_FIELD_SPEC.md): the discrete specification of the compact correction-field (CCF) coupling.

The default nested coupling is PTF (`crescent.l1_coupling_mode=1`); CCF is
mode 2.

## Sample Cases

- [test_case/README.md](../test_case/README.md): curated sample input cases.

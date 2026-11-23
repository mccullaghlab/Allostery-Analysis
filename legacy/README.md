# Original code (reference only)

This directory contains the code as used for the 2020 JCTC paper. It is **not** used by the
`allostery` package. It is kept as a reference and for the regression tests in
`tests/test_legacy_regression.py`, which compile these programs for a small system and
check that the new implementation reproduces their output exactly.

* `henm/`: hENM fitting scripts. Note that `henm_functions.py` contains a bug: the
  distance-cutoff mask is inverted. It is fixed in `allostery.henm`.
* `paths/`: Fortran path sampling (`paths_TWISP*.f`), analysis programs
  (`paths_analyze*`, `paths_distill*`, `paths_deg*`, `paths_partition*`), the driver
  script `run_multi.sh`, and the derivative scripts `edge_effect*_multi.py`. All have the
  system size (`npar=454`) and file names hard-coded.

The pre-compiled macOS binaries (`*.e`) that used to be committed have been removed.
Build the programs with `make` in `paths/` if needed.

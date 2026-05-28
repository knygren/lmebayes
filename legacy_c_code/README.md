# Legacy vendored R Mathlib (CPU)

Archive of **src/** C sources and headers that duplicated R’s Mathlib
(`src/nmath` subset compiled into the package DLL). C++ now calls libR via
`<Rmath.h>` (`Rf_dnorm4`, `Rf_pgamma`, `Rf_qgamma`, `Rf_dbinom_raw`, etc.) in:

- `famfuncs_gaussian.cpp`, `famfuncs_Gamma.cpp`, `famfuncs_binomial.cpp`
- `rng_utils.cpp`

## Contents

- **\*.c** — former `src/*.c` translation units (bd0, dbinom, dnorm, …)
- **\*.h** — local stand-ins for `Rmath.h` / `dpq.h` / `Rconfig` used by those `.c` files

## Not archived here

- **`src/nmath/`** — full Mathlib snapshot (reference / OpenCL port source); not linked into the R package CPU build
- **`inst/cl/nmath/`** — OpenCL ports (still required for GPU)

Safe to remove the archived files from **`src/`** once the package builds and tests pass without them.

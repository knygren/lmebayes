# `R/` — exported and documented functions

Symbols defined in **`R/`** that are **exported** (`NAMESPACE`) or have a **help
page** (`man/*.Rd`). Use this list when reviewing the public API, `\usage`
blocks, and README coverage.

Companion: [R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md) (`@noRd` and other
undocumented helpers).

User-facing export tables also appear in the [README Function overview](../README.md#function-overview).

Mixed-model setup (`model_setup`, `Prior_Setup_lmebayes`, lme4 design utilities,
`rlmerb` / `rglmerb` engines) live in **glmbayesCore**; see
[glmbayesCore/inst/R_EXPORTED_AND_DOCUMENTED.md](../../glmbayesCore/inst/R_EXPORTED_AND_DOCUMENTED.md).

---

## Native exports (`R/` → `NAMESPACE`)

| Function | File | Role |
|----------|------|------|
| `lmerb()` | `lmerb.R` | Formula LMM; calls `rlmerb()`. |
| `glmerb()` | `glmerb.R` | Formula GLMM; calls `rglmerb()`. |
| `print.lmerb()` | `lmerb.R` | Fit print; optional sweep-history tables. |
| `print.glmerb()` | `glmerb.R` | Same for GLMM. |
| `summary.lmerb()` | `summary.lmerb.R` | Detailed posterior / prior / ranef summary. |
| `summary.glmerb()` | `summary.lmerb.R` | Alias of `summary.lmerb`. |
| `print.summary.lmerb()` | `summary.lmerb.R` | Print method for summary object. |
| `Prior_SetupBlock()` | `prior_setup_lmebayes.R` | Per-row-block prior setup. |
| `lmbBlock()` | `lmbBlock.R` | Row-block Bayesian LM fits. |
| `glmbBlock()` | `glmbBlock.R` | Row-block Bayesian GLM fits. |
| `print.blmb()` | `lmbBlock.R` | Print block LM list. |
| `print.bglmb()` | `glmbBlock.R` | Print block GLM list. |
| `summary.blmb()` / `print.summary.blmb()` | `summary.blmb.R` | Block LM summaries. |
| `summary.bglmb()` / `print.summary.bglmb()` | `summary.bglmb.R` | Block GLM summaries. |
| `has_opencl()` | `has_opencl.R` | Compile-time OpenCL flag for this build. |

---

## Re-exports (`R/reexports_*.R`)

Implemented in **glmbayes** / **glmbayesCore**; documented under dependency help pages.

| Function | Source package |
|----------|----------------|
| `lmb()`, `glmb()`, `directional_tail()` | glmbayes |
| `Prior_Setup()`, `dNormal()`, `dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, `dGamma()` | glmbayesCore |
| `pfamily_list()` | glmbayesCore |
| `plot_sweep_history_diag()` | glmbayesCore |
| `model_setup()`, `Prior_Setup_lmebayes()`, `rlmerb()`, `rglmerb()` | glmbayesCore |

S3 `print.model_setup` and `print.lmebayes_prior_setup` register in **glmbayesCore**;
**lmebayes** dispatches them via `import(glmbayesCore)`.

---

## Documented but not exported (`@keywords internal` + `man/`)

Callable with `lmebayes:::`; have help pages but are not in `NAMESPACE`.

| Function | File | Role | Review |
|----------|------|------|--------|
| `print_coef_means()` | `glmerb.R` | MLE vs mode vs posterior mean table. | **Export?** Heavily used in demos. |
| `block_check_identifiability()` | `lmbBlock.R` | Block GLM rank (formula). | Block models. |
| `block_check_identifiability_xy()` | `lmbBlock.R` | Block GLM rank (matrix). | Block models. |

lme4 design utilities (`get_lme4_components`, `extract_re_hyper_matrices`, …) moved to
**glmbayesCore** (`R/lme4_design_utilities.R`); see Core `inst/R_EXPORTED_AND_DOCUMENTED.md`.

---

## Documentation topics (no function body in `R/`)

| Topic / file | Contents |
|--------------|----------|
| `lmebayes-package.R` | Package meta, imports (`"_PACKAGE"`). |
| `gpu_diagnostics.R` | `@name gpu_diagnostics` — links to OpenCL diagnostics. |
| `data-*.R` | Lazy data docs (`Boston_centered`, `BikeSharing`, `carinsca`, etc.). |

---

## Gaps (not in `R/` yet)

Ex. 16 demo logic (proxy \(\hat D_\ell\), variance ratios) is **not** exported or
documented here. Candidate home: extend `plot_sweep_history_diag()` or add
`sweep_history_diag_*()` in `R/`.

---

## Review checklist (exports / docs)

| Priority | Item |
|----------|------|
| 1 | Export `print_coef_means()` or document `:::` as intentional in README. |
| 2 | Extend sweep-history diagnostics in `R/` when Ex. 16 probes stabilize. |
| 3 | Run `devtools::document()` after any `@export` or `\usage` change. |

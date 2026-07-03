# `R/` — exported and documented functions

Symbols defined in **`R/`** that are **exported** (`NAMESPACE`) or have a **help
page** (`man/*.Rd`). Use this list when reviewing the public API, `\usage`
blocks, and README coverage.

Companion: [R_INTERNAL_HELPERS.md](R_INTERNAL_HELPERS.md) (`@noRd` and other
undocumented helpers).

User-facing export tables also appear in the [README Function overview](../README.md#function-overview).

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
| `rlmerb()` | `rlmerb.R` | Matrix-level Gaussian LMM sampler. |
| `rglmerb()` | `rglmerb.R` | Matrix-level GLMM sampler. |
| `model_setup()` | `model_setup.R` | Parse lme4-style formula → design object. |
| `print.model_setup()` | `model_setup.R` | Print design object. |
| `Prior_Setup_lmebayes()` | `prior_setup_lmebayes.R` | Calibrate Block 2 priors from lmer/glmer. |
| `Prior_SetupBlock()` | `prior_setup_lmebayes.R` | Per-row-block prior setup. |
| `print.lmebayes_prior_setup()` | `prior_setup_lmebayes.R` | Print prior-setup object. |
| `pfamily_list()` method | `pfamily_list.R` | S3 for `lmebayes_prior_setup`. |
| `lmbBlock()` | `lmbBlock.R` | Row-block Bayesian LM fits. |
| `glmbBlock()` | `glmbBlock.R` | Row-block Bayesian GLM fits. |
| `print.blmb()` | `lmbBlock.R` | Print block LM list. |
| `print.bglmb()` | `glmbBlock.R` | Print block GLM list. |
| `summary.blmb()` / `print.summary.blmb()` | `summary.blmb.R` | Block LM summaries. |
| `summary.bglmb()` / `print.summary.bglmb()` | `summary.bglmb.R` | Block GLM summaries. |
| `plot_sweep_history_diag()` | `plot_sweep_history_diag.R` | Cross-chain mean/SD vs inner sweep. |
| `has_opencl()` | `has_opencl.R` | Compile-time OpenCL flag for this build. |

---

## Re-exports (`R/reexports_*.R`)

Implemented in **glmbayes** / **glmbayesCore**; documented under `?reexports` and
dependency help pages.

| Function | Source package |
|----------|----------------|
| `lmb()`, `glmb()`, `directional_tail()` | glmbayes |
| `Prior_Setup()`, `dNormal()`, `dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, `dGamma()` | glmbayesCore |

---

## Documented but not exported (`@keywords internal` + `man/`)

Callable with `lmebayes:::`; have help pages but are not in `NAMESPACE`.

| Function | File | Role | Review |
|----------|------|------|--------|
| `print_coef_means()` | `glmerb.R` | MLE vs mode vs posterior mean table. | **Export?** Heavily used in demos. |
| `is_single_factor_model()` | `glmerb_utilities.R` | Exactly one grouping factor. | Keep internal. |
| `is_fixed_effects_only()` | `glmerb_utilities.R` | No random effects in formula. | Keep internal. |
| `get_lme4_components()` | `glmerb_utilities.R` | lme4 parse → matrices. | Dev / extension utility. |
| `show_lme4_Z_random()` | `glmerb_utilities.R` | Debug random-effects design. | Dev only. |
| `classify_lme4_fixed_columns()` | `glmerb_utilities.R` | Population vs group-level fixed cols. | Prior pipeline. |
| `classify_crosslevel_re_moderation()` | `glmerb_utilities.R` | Cross-level RE structure. | Prior pipeline. |
| `extract_re_hyper_matrices()` | `glmerb_utilities.R` | Block 2 group designs. | Prior / setup. |
| `extract_re_Z_obs()` | `glmerb_utilities.R` | Obs-level Z for one group. | Prior / setup. |
| `extract_lme4_submatrices()` | `glmerb_utilities.R` | Subset parsed lme4 parts. | Prior / setup. |
| `extract_lme4_fixed_group_matrix()` | `glmerb_utilities.R` | Group-level fixed matrix. | Prior / setup. |
| `extract_lmer_variance_components()` | `glmerb_utilities.R` | From **lmer** fit. | Prior calibration. |
| `extract_mer_variance_components()` | `glmerb_utilities.R` | From **glmer** fit. | Prior calibration. |
| `lmerb_default_vcov_formula()` | `glmerb_utilities.R` | Default vcov for prior scaling. | Prior calibration. |
| `block_check_identifiability()` | `lmbBlock.R` | Block GLM rank (formula). | Block models. |
| `block_check_identifiability_xy()` | `lmbBlock.R` | Block GLM rank (matrix). | Block models. |

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

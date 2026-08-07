# `R/` — internal helpers (undocumented in `man/`)

Functions and symbols in **`R/`** with **`@noRd`**, no roxygen, or
`@keywords internal` **without** a dedicated help page. Intended for
`lmebayes:::` only — not part of the exported API.

**Columns:** *File* is the defining source; *Called from* lists direct callers
in `R/` (comma-separated). Helpers with no callers are marked *(unused)*.

Companion: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).

Mixed-model sampling, prior-setup, and lme4 design helpers were moved to
**glmbayesCore**; see
[glmbayesCore/inst/R_INTERNAL_HELPERS.md](../../glmbayesCore/inst/R_INTERNAL_HELPERS.md).

---

## Mixed-model glue (`glmerb_utilities.R`)

| Function | File | Called from |
|----------|------|-------------|
| `.lmerb_reference_fit()` | `glmerb_utilities.R` | `summary.lmerb()`, `.lmerb_fixef_component_summary()`, `.lmerb_tau2_prior_overview()` |
| `.lmebayes_stage_v2_fixef()` | `glmerb_utilities.R` | *(unused; legacy two-block staging)* |

---

## Resolved from **glmbayesCore** (`glmbayesCore:::` / `importFrom`)

| Symbol | Core file | Called from (lmebayes) |
|--------|-----------|------------------------|
| `.lmebayes_priors_from_pfamily_list()` | `mixed_rmerb_helpers.R` | `lmerb()`, `glmerb()` |
| `.lmebayes_block2_icm_labels()` | `mixed_rmerb_helpers.R` | `lmerb()`, `glmerb()` |
| `.lmebayes_mer_optional_args()` | `model_setup.R` | `glmerb()` |
| `extract_mer_variance_components()` | `lme4_design_utilities.R` | `summary.lmerb()` |
| `.two_block_as_staged_names()` | `two_block_glmm_pilot_helpers.R` | `.lmebayes_stage_v2_fixef()` |
| `.mrglmb_normalize_pfamily_lists` | *(Core)* | `lmbBlock()`, `glmbBlock()` (alias in `block_core_pfamily.R`) |
| `.validate_pfamily_for_rlmb` | *(Core)* | `lmbBlock()`, `glmbBlock()` (alias in `block_core_pfamily.R`) |

Formula drivers call re-exported `model_setup()`, `Prior_Setup_GLMM()`, `rlmerb()`,
and `rglmerb()` without namespace qualification.

**Direct Core calls (must stay exported):** `build_mu_all()`, `lmerb_posterior_mean()`,
`glmerb_posterior_mode()` (`importFrom`). See
[R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).

**Direct `lmebayesCore` call:** `normalize_block()` (`lmebayesCore::` in
`.blmb_formula_block_meta()`; also called directly by `block_check_identifiability_xy()`).
This is `lmbBlock()`/`glmbBlock()`/`Prior_SetupBlock()`'s only runtime dependency on
`lmebayesCore` — full trace in
[lmebayesCore/inst/LMBBLOCK_LMEBAYESCORE_DEPENDENCIES.md](../../lmebayesCore/inst/LMBBLOCK_LMEBAYESCORE_DEPENDENCIES.md).

**Indirect only (export optional for lmebayes):** `rGLMM()`, `rLMMNormal_reg()`,
`rLMMNormal_reg_estimated_vcov()`, `rLMMindepNormalGamma_reg()` — listed under
**glmbayesCore-only exports** in Core inventory.

---

## `summary.lmerb()` helpers (`summary.lmerb.R`)

| Function | File | Called from |
|----------|------|-------------|
| `.lmerb_print_summary_table()` | `summary.lmerb.R` | `print.summary.lmerb()` |
| `.lmerb_lmer_fixef_lookup()` | `summary.lmerb.R` | `.lmerb_fixef_component_summary()` |
| `.lmerb_fixef_component_summary()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_fixef_prior_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_fixef_percentiles_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_fixef_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_tau2_prior_overview()` | `summary.lmerb.R` | `summary.lmerb()`, `.lmerb_tau2_posterior_overview()`, `.lmerb_tau2_sd_percentiles_overview()` |
| `.lmerb_tau2_posterior_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_tau2_percentiles_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_tau2_sd_percentiles_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_ranef_overview()` | `summary.lmerb.R` | `summary.lmerb()` |
| `.lmerb_ranef_groups_detail()` | `summary.lmerb.R` | `summary.lmerb()` |

---

## Row-block helpers (`lmbBlock.R`, `glmbBlock.R`)

| Function | File | Called from |
|----------|------|-------------|
| `.blmb_formula_block_meta()` | `lmbBlock.R` | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()`, `.blmb_blocks_full_rank()` |
| `.blmb_rows_to_data_subset()` | `lmbBlock.R` | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()` |
| `.blmb_resolve_block()` | `lmbBlock.R` | `.blmb_formula_block_meta()` |
| `.blmb_lmb_display_call()` | `lmbBlock.R` | `lmbBlock()` |
| `.blmb_assemble()` | `lmbBlock.R` | `lmbBlock()` |
| `.blmb_blocks_full_rank()` | `lmbBlock.R` | `block_check_identifiability()` |
| `.blmb_blocks_full_rank_xy()` | `lmbBlock.R` | `block_check_identifiability_xy()` |
| `.blmb_coef_means_matrix()` | `lmbBlock.R` | `summary.blmb()`, `summary.bglmb()`, `print.blmb()`, `print.bglmb()` |
| `.blmb_dic_table()` | `lmbBlock.R` | `summary.blmb()`, `summary.bglmb()`, `print.blmb()`, `print.bglmb()` |
| `.blmb_glmb_display_call()` | `glmbBlock.R` | `glmbBlock()` |
| `.bglmb_assemble()` | `glmbBlock.R` | `glmbBlock()` |
| `.blmb_resolve_block_calibration_arg()` | `prior_setup_lmebayes.R` | `Prior_SetupBlock()` |

Full inventory and call graphs, including the `lmbBlock()`/`glmbBlock()`
symmetry and where `Prior_SetupBlock()` fits in:
[README_LMBBLOCK.md](README_LMBBLOCK.md).

---

## Build, attach, C++

| Function | File | Called from |
|----------|------|-------------|
| `.has_opencl_cpp()` | `rcpp_wrappers.R` | `has_opencl()` |
| `has_opencl_cpp_export()` | `RcppExports.R` | `.has_opencl_cpp()` (`.Call` entry) |
| `use_RcppParallel()` | `internal_rcppparallel.R` | *(unused at runtime; satisfies R CMD check `Imports: RcppParallel`)* |
| `.opencl_startup_quiet()` | `zzz.R` | `.opencl_startup_message()` |
| `.opencl_runtime_sniff()` | `zzz.R` | `.opencl_startup_message()` |
| `.opencl_startup_message()` | `zzz.R` | `.onAttach()` |
| `.onAttach()` | `zzz.R` | R load hook (calls `.opencl_startup_message()`) |

---

## Load-time aliases (`block_core_pfamily.R`)

| Symbol | File | Called from |
|--------|------|-------------|
| `.mrglmb_normalize_pfamily_lists` | `block_core_pfamily.R` | `lmbBlock()`, `glmbBlock()` |
| `.validate_pfamily_for_rlmb` | `block_core_pfamily.R` | `lmbBlock()`, `glmbBlock()` (passed to `.mrglmb_normalize_pfamily_lists()`) |

Both symbols are namespace aliases to **glmbayesCore** helpers.

---

## Review checklist (internals)

| Priority | Item |
|----------|------|
| 1 | Avoid adding new `@noRd` helpers unless tied to exported behavior. |
| 2 | Promote repeated demo patterns to documented `R/` API instead of new `:::` helpers. |
| 3 | Remove or wire up *(unused)* helpers (`.lmebayes_stage_v2_fixef`) when touching related code. |

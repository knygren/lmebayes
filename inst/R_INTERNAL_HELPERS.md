# `R/` — internal helpers (undocumented in `man/`)

Functions and symbols in **`R/`** with **`@noRd`**, no roxygen, or
`@keywords internal` **without** a dedicated help page. Intended for
`lmebayes:::` only — not part of the exported API.

**Columns:** *File* is the defining source; *Called from* lists direct callers
in `R/` (comma-separated). Helpers with no callers are marked *(unused)*.

Companion: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).

---

## Sampling engines (`glmerb_utilities.R`)

| Function | File | Called from |
|----------|------|-------------|
| `.lmebayes_run_lmm_engine()` | `glmerb_utilities.R` | `rglmerb()` (Gaussian path), `rlmerb()` |
| `.lmebayes_add_fixef_summaries()` | `glmerb_utilities.R` | `rglmerb()`, `rlmerb()` |
| `.lmebayes_stage_v2_fixef()` | `glmerb_utilities.R` | *(unused; legacy two-block staging)* |

---

## Priors / dispersion

| Function | File | Called from |
|----------|------|-------------|
| `.lmebayes_priors_from_pfamily_list()` | `glmerb_utilities.R` | `lmerb()`, `glmerb()` |
| `.lmebayes_resolve_dispersion_ranef()` | `glmerb_utilities.R` | `rglmerb()`, `rlmerb()`, `.lmebayes_priors_from_pfamily_list()` |
| `.lmebayes_validate_dispersion_ranef()` | `glmerb_utilities.R` | *(unused; thin wrapper around `.lmebayes_resolve_dispersion_ranef()`)* |
| `.lmebayes_block1_prior_list()` | `glmerb_utilities.R` | `rglmerb()`, `rlmerb()` |
| `.lmerb_reference_fit()` | `glmerb_utilities.R` | `summary.lmerb()`, `.lmerb_fixef_component_summary()`, `.lmerb_tau2_prior_overview()` |
| `.lmebayes_resolve_pwt()` | `prior_setup_lmebayes.R` | `Prior_Setup_lmebayes()` |
| `.lmebayes_resolve_disp_prior()` | `prior_setup_lmebayes.R` | `Prior_Setup_lmebayes()` |
| `.lmebayes_block_glm_estimable()` | `prior_setup_lmebayes.R` | `Prior_Setup_lmebayes()` |

---

## Formula / design parsing

| Function | File | Called from |
|----------|------|-------------|
| `.lmebayes_validate_uncorrelated_re_formula()` | `glmerb_utilities.R` | `extract_re_hyper_matrices()` |
| `.lmebayes_mer_convergence_issues()` | `model_setup.R` | `Prior_Setup_lmebayes()` |
| `.lmebayes_normalize_family()` | `model_setup.R` | `model_setup()` |
| `.lmebayes_mer_optional_args()` | `model_setup.R` | `model_setup()`, `glmerb()` |
| `.lme4_Z_random_column_map()` | `glmerb_utilities.R` | `.lme4_Z_random_colnames()`, `get_lme4_components()` |
| `.lme4_Z_random_colnames()` | `glmerb_utilities.R` | `.lme4_label_Z_random_sparse()` |
| `.lme4_Z_random_rownames()` | `glmerb_utilities.R` | `.lme4_Z_random_row_map()`, `.lme4_label_Z_random_sparse()` |
| `.lme4_Z_random_row_map()` | `glmerb_utilities.R` | `get_lme4_components()` |
| `.lme4_label_Z_random_sparse()` | `glmerb_utilities.R` | `get_lme4_components()` |

Exported entry points that reach the Z-label chain: `model_setup()` →
`extract_re_hyper_matrices()` → `get_lme4_components()`.

---

## Sampling console output (`rglmerb_diag.R`)

| Function | File | Called from |
|----------|------|-------------|
| `.lmebayes_block2_icm_labels()` | `rglmerb_diag.R` | `rglmerb()`, `rlmerb()`, `lmerb()` (ICM-only path), `glmerb()` (ICM-only path) |
| `.lmebayes_print_icm_fixef_table()` | `rglmerb_diag.R` | `rglmerb()`, `rlmerb()` |
| `.lmebayes_print_ranef_mode_reference()` | `rglmerb_diag.R` | `rglmerb()` |
| `.lmebayes_print_fixef_init()` | `rglmerb_diag.R` | `rglmerb()`, `rlmerb()` |

`lmerb()` and `glmerb()` inline the ICM fixed-effect table using
`.lmebayes_block2_icm_labels()` only (no `print_icm_fixef_table()`).

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
| 3 | Remove or wire up *(unused)* helpers (`.lmebayes_stage_v2_fixef`, `.lmebayes_validate_dispersion_ranef`) when touching related code. |

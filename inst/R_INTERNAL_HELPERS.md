# `R/` — internal helpers (undocumented in `man/`)

Functions and symbols in **`R/`** with **`@noRd`**, no roxygen, or
`@keywords internal` **without** a dedicated help page. Intended for
`lmebayes:::` only — not part of the exported API.

Companion: [R_EXPORTED_AND_DOCUMENTED.md](R_EXPORTED_AND_DOCUMENTED.md).

---

## Sampling engines

| Function | File | Called from |
|----------|------|-------------|
| `.rglmerb_v6_rGLMM()` | `rglmerb_v6.R` | `rglmerb()` non-Gaussian path |
| `.lmebayes_run_lmm_engine()` | `glmerb_utilities.R` | `rglmerb()` Gaussian path |
| `.lmebayes_stage_v2_fixef()` | `glmerb_utilities.R` | Engine staging |
| `.lmebayes_add_fixef_summaries()` | `glmerb_utilities.R` | Post-process sampler output |

**Rename note:** `rglmerb_v6.R` / `.rglmerb_v6_rGLMM()` — optional future rename
to drop the `v6` suffix (legacy `rglmerb_v5` source removed).

---

## Priors / dispersion

| Function | File |
|----------|------|
| `.lmebayes_priors_from_pfamily_list()` | `glmerb_utilities.R` |
| `.lmebayes_resolve_pwt()` | `prior_setup_lmebayes.R` |
| `.lmebayes_resolve_disp_prior()` | `prior_setup_lmebayes.R` |
| `.lmebayes_resolve_dispersion_ranef()` | `glmerb_utilities.R` |
| `.lmebayes_validate_dispersion_ranef()` | `glmerb_utilities.R` |
| `.lmebayes_block_glm_estimable()` | `prior_setup_lmebayes.R` |
| `.lmebayes_block1_prior_list()` | `glmerb_utilities.R` |
| `.lmerb_reference_fit()` | `glmerb_utilities.R` |

---

## Formula / design parsing

| Function | File |
|----------|------|
| `.lmebayes_validate_uncorrelated_re_formula()` | `glmerb_utilities.R` |
| `.lmebayes_mer_convergence_issues()` | `model_setup.R` |
| `.lmebayes_normalize_family()` | `model_setup.R` |
| `.lmebayes_mer_optional_args()` | `model_setup.R` |
| `.lme4_Z_random_column_map()` | `glmerb_utilities.R` |
| `.lme4_Z_random_colnames()` | `glmerb_utilities.R` |
| `.lme4_Z_random_rownames()` | `glmerb_utilities.R` |
| `.lme4_Z_random_row_map()` | `glmerb_utilities.R` |
| `.lme4_label_Z_random_sparse()` | `glmerb_utilities.R` |

---

## Sampling console output

| Function | File |
|----------|------|
| `.lmebayes_block2_icm_labels()` | `rglmerb_diag.R` |
| `.lmebayes_print_icm_fixef_table()` | `rglmerb_diag.R` |
| `.lmebayes_print_ranef_mode_reference()` | `rglmerb_diag.R` |
| `.lmebayes_print_fixef_init()` | `rglmerb_diag.R` |

---

## `summary.lmerb()` helpers (`summary.lmerb.R`)

| Function | Role |
|----------|------|
| `.lmerb_print_summary_table()` | Format summary tables. |
| `.lmerb_lmer_fixef_lookup()` | Map RE/param to lmer fixed effects. |
| `.lmerb_fixef_component_summary()` | One RE component block. |
| `.lmerb_fixef_prior_overview()` | Prior slice of fixef summary. |
| `.lmerb_fixef_percentiles_overview()` | Fixef percentiles. |
| `.lmerb_fixef_overview()` | Top-level fixef section. |
| `.lmerb_tau2_prior_overview()` | ING τ² prior section. |
| `.lmerb_tau2_posterior_overview()` | ING τ² posterior section. |
| `.lmerb_tau2_percentiles_overview()` | τ² percentiles. |
| `.lmerb_tau2_sd_percentiles_overview()` | τ² SD percentiles. |
| `.lmerb_ranef_overview()` | Ranef summary section. |
| `.lmerb_ranef_groups_detail()` | Per-group ranef detail. |

---

## Row-block helpers (`lmbBlock.R`, `glmbBlock.R`)

| Function | File | Role |
|----------|------|------|
| `.blmb_formula_block_meta()` | `lmbBlock.R` | Block metadata from formula. |
| `.blmb_rows_to_data_subset()` | `lmbBlock.R` | Row index → data subset. |
| `.blmb_resolve_block()` | `lmbBlock.R` | Resolve block factor. |
| `.blmb_lmb_display_call()` | `lmbBlock.R` | Display call for one block. |
| `.blmb_assemble()` | `lmbBlock.R` | Assemble block fit list. |
| `.blmb_coef_means_matrix()` | `lmbBlock.R` | Coef means across blocks. |
| `.blmb_dic_table()` | `lmbBlock.R` | DIC table for blocks. |
| `.blmb_blocks_full_rank()` | `lmbBlock.R` | Full-rank check. |
| `.blmb_blocks_full_rank_xy()` | `lmbBlock.R` | Full-rank check (matrix). |
| `.blmb_glmb_display_call()` | `glmbBlock.R` | GLM block display call. |
| `.bglmb_assemble()` | `glmbBlock.R` | Assemble GLM block list. |

---

## Build, attach, C++

| Function | File |
|----------|------|
| `.has_opencl_cpp()` | `rcpp_wrappers.R` |
| `has_opencl_cpp_export()` | `RcppExports.R` |
| `use_RcppParallel()` | `internal_rcppparallel.R` |
| `.opencl_startup_quiet()` | `zzz.R` |
| `.opencl_runtime_sniff()` | `zzz.R` |
| `.opencl_startup_message()` | `zzz.R` |
| `.onAttach()` | `zzz.R` |

---

## Load-time aliases (`block_core_pfamily.R`)

| Symbol | Role |
|--------|------|
| `.mrglmb_normalize_pfamily_lists` | Alias to **glmbayesCore** |
| `.validate_pfamily_for_rlmb` | Alias to **glmbayesCore** |

---

## Review checklist (internals)

| Priority | Item |
|----------|------|
| 1 | Rename `rglmerb_v6.R` when version suffix is dropped. |
| 2 | Avoid adding new `@noRd` helpers unless tied to exported behavior. |
| 3 | Promote repeated demo patterns to documented `R/` API instead of new `:::` helpers. |

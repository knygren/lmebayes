# Row-block models — full function reference

Everything behind the three row-block entry points (SAS `BY`-style splits:
independent fits/priors per group, sharing one formula). Expands the
`README.md` "Row-block models" / "Model setup and mixed-model priors"
tables into a complete inventory, and explains the `lmbBlock()`/`glmbBlock()`
symmetry underneath them.

Supersedes an earlier `lmbBlock()`-centric draft of this file with a
symmetric treatment of all three entry points.

**Last reviewed:** 2026-07-20.

---

## At a glance

### Row-block models (SAS `BY`-style splits)

| Function | Role |
|----------|------|
| `lmbBlock()` | One `glmbayes::lmb()` fit per row block (shared formula, block-specific priors). |
| `glmbBlock()` | One `glmbayes::glmb()` fit per row block. |
| `print()` / `summary()` | S3 methods for `"blmb"` and `"bglmb"` (lists of block fits). |

### Model setup and mixed-model priors

| Function | Role |
|----------|------|
| `Prior_SetupBlock()` | Run `glmbayesCore::Prior_Setup()` independently on each row block (for `lmbBlock()` / `glmbBlock()`). |

All three top-level functions — `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()`
— follow the exact same shape: *reconcile the `block` argument against
`model.frame()`, loop over blocks calling one external single-fit function
per block, then assemble the `k` results into one object.* They differ only
in which single-fit function they loop (`lmb()`, `glmb()`, or `Prior_Setup()`)
and what the assembled result needs to carry.

---

## The `lmbBlock()` / `glmbBlock()` symmetry

`glmbBlock()` is `lmbBlock()`'s GLM counterpart, built the same way `glmb()`
is `lmb()`'s GLM counterpart in `glmbayes` itself. Concretely:

### Literally shared (same function object, called by both)

| Function | Role |
|----------|------|
| `.blmb_formula_block_meta()` | `model.frame()`/`model.matrix()` once, reconcile `block` against it, call `normalize_block()`. |
| `.blmb_resolve_block()` | Accepts the 5 `block` input shapes (factor / vector / column name / list / `l2_blocks`). |
| `.blmb_rows_to_data_subset()` | Translates model-frame row indices back to original-`data` row indices/rownames for `subset=`. |
| `.mrglmb_normalize_pfamily_lists()` (from `glmbayesCore`) | Validate/recycle `pfamily` vs. `pfamily_list` into one-per-block. |
| `.validate_pfamily_for_rlmb()` (from `glmbayesCore`) | Per-block pfamily validator passed into the above. |
| `.blmb_coef_means_matrix()` | Stack per-block posterior mean coefficients into one matrix (used by both `print`/`summary` pairs). |
| `.blmb_dic_table()` | Stack per-block `DIC`/`pD` into one table (used by both `print`/`summary` pairs). |
| `lmebayesCore::normalize_block()` | Canonicalize the reconciled `block` vector into `list(k, ids, l2_blocks, starts, rows)`. |

### Mirrored (parallel, near-identical pairs — one per family)

| `lmbBlock` side | `glmbBlock` side | What differs |
|-----------------|--------------------|----------------|
| `glmbayes::lmb()` | `glmbayes::glmb()` | Different engine — the actual statistical difference; everything else below exists only to carry this choice through. |
| `.blmb_lmb_display_call()` | `.blmb_glmb_display_call()` | Builds `call("lmb", ...)` vs `call("glmb", ...)`; the `glmb` version's pass-through arg list adds `"family"`. |
| `.blmb_assemble()` | `.bglmb_assemble()` | `bglmb` version additionally stores `attr(outlist, "family")`. |
| `print.blmb()` | `print.bglmb()` | `bglmb` version has an extra block printing `family$family`/`family$link`; otherwise identical, and both call the *same* `.blmb_coef_means_matrix()` / `.blmb_dic_table()`. |
| `summary.blmb()` / `print.summary.blmb()` | `summary.bglmb()` / `print.summary.bglmb()` | Same pattern: `bglmb` versions carry/print `family`, otherwise identical. |

There's also a small asymmetry in the *public* argument lists (not a backend
function difference): `lmbBlock()` exposes `method = "qr"`, `qr`,
`singular.ok` (`lm()`-only concepts), while `glmbBlock()` exposes `family`
instead — reflecting `lm()` vs `glm()`, not a design inconsistency.

**Why mirror instead of parametrize a single implementation?** The
alternative would be one function branching on `is.null(family)` throughout
— `.blmb_lmb_display_call()`/`.blmb_glmb_display_call()`,
`.blmb_assemble()`/`.bglmb_assemble()`, and both `print`/`summary` pairs would
each become one function with conditional logic instead of two simple ones.
`glmbayes` itself made the same choice (`lmb()` and `glmb()` are separate
functions, not one function branching on family), so `lmbBlock()`/`glmbBlock()`
mirror that precedent rather than diverging from it one layer up.

---

## `Prior_SetupBlock()`: the third sibling

`Prior_SetupBlock()` reuses the **same shared plumbing** as `lmbBlock()`/`glmbBlock()`
(`.blmb_formula_block_meta()`, `.blmb_rows_to_data_subset()`) — this is the
reason those two helpers are standalone functions rather than being inlined
into `lmbBlock()`. It loops `glmbayesCore::Prior_Setup()` (re-exported by
`lmebayes`) once per block instead of `lmb()`/`glmb()`, and has one helper of
its own that neither `lmbBlock()` nor `glmbBlock()` need:

| Function | Role |
|----------|------|
| `.blmb_resolve_block_calibration_arg()` | Lets each of `Prior_Setup()`'s six calibration args (`pwt`, `n_prior`, `sd`, `mu`, `dispersion`, `k`) be supplied either as one shared value (recycled to every block, `Prior_Setup()`'s own default behavior) or as a named list keyed by block ID (one value per block, validated for exact name match against `block_ids`). |

This exists because `Prior_Setup()`'s calibration knobs are genuinely
useful to vary *per block* (e.g. a tighter prior `sd` for a block with less
data) — something `lmb()`/`glmb()`'s own arguments don't need an analogous
per-block override for, since block-specific priors there are already
expressed via `pfamily_list` (one full `pfamily` object per block, not a
scalar knob).

---

## Optional preflight diagnostics (not called by any of the three)

Documented alongside `lmbBlock()` but standalone, not part of any of the
three call graphs above — users invoke these separately, before fitting:

| Function | Exported | Role |
|----------|----------|------|
| `block_check_identifiability()` | No (`lmebayes:::`) | Two-level rank preflight, formula interface: Level 1 = each block's design matrix full column rank; Level 2 = hyper design (`X_nbhd`) full rank across Level-1-identified blocks. |
| `block_check_identifiability_xy()` | No (`lmebayes:::`) | Same, matrix interface (for use inside a Gibbs-loop setup); also calls `lmebayesCore::normalize_block()` directly. |
| `.blmb_blocks_full_rank()` | No | Level-1 rank table, formula interface. |
| `.blmb_blocks_full_rank_xy()` | No | Level-1 rank table, matrix interface. |

For plain `lmbBlock()`/`glmbBlock()` (independent per-block fits, no shared
population parameter), only Level 1 is meaningful — Level 2 is really for
the *coupled* Gibbs mixed-model samplers (`lmerb()`/`glmerb()`), where
identifying a shared `mu` across blocks genuinely matters. The same
two-level algorithm is exposed here anyway since it costs nothing extra and
keeps one implementation instead of two.

---

## Call graphs

```
lmbBlock()
├── .blmb_formula_block_meta()
│   ├── stats::model.frame() / model.matrix()
│   ├── .blmb_resolve_block()
│   └── lmebayesCore::normalize_block()
├── .mrglmb_normalize_pfamily_lists()  [from glmbayesCore]
│   └── .validate_pfamily_for_rlmb()   [from glmbayesCore]
├── (loop b = 1..k)
│   ├── .blmb_rows_to_data_subset()
│   ├── do.call(glmbayes::lmb, ...)
│   └── .blmb_lmb_display_call()
└── .blmb_assemble()

glmbBlock()
├── .blmb_formula_block_meta()          [same as lmbBlock()]
├── .mrglmb_normalize_pfamily_lists()   [same as lmbBlock()]
├── (loop b = 1..k)
│   ├── .blmb_rows_to_data_subset()     [same as lmbBlock()]
│   ├── do.call(glmbayes::glmb, ...)    ← the actual difference
│   └── .blmb_glmb_display_call()       ← mirrored, not shared
└── .bglmb_assemble()                   ← mirrored, not shared

Prior_SetupBlock()
├── .blmb_formula_block_meta()          [same as lmbBlock()/glmbBlock()]
├── .blmb_resolve_block_calibration_arg()  [× 6 calibration args; unique to this entry point]
└── (loop b = 1..k)
    ├── .blmb_rows_to_data_subset()     [same as lmbBlock()/glmbBlock()]
    └── do.call(glmbayesCore::Prior_Setup, ...)

print()/summary() on a "blmb" result
└── .blmb_coef_means_matrix(), .blmb_dic_table()   [shared with "bglmb"]

print()/summary() on a "bglmb" result
└── .blmb_coef_means_matrix(), .blmb_dic_table()   [shared with "blmb"]
    + family/link printing                          [mirrored addition]
```

---

## Full function inventory

| Function | File | Exported | Shared with | Role |
|----------|------|----------|-------------|------|
| `lmbBlock()` | `lmbBlock.R` | Yes | — | Main driver, Gaussian (`lmb()` per block). |
| `glmbBlock()` | `glmbBlock.R` | Yes | — | Main driver, any GLM family (`glmb()` per block). |
| `Prior_SetupBlock()` | `prior_setup_lmebayes.R` | Yes | — | `Prior_Setup()` per block. |
| `print.blmb()` | `lmbBlock.R` | Yes (S3) | mirrors `print.bglmb()` | Combined print for `lmbBlock()` result. |
| `print.bglmb()` | `glmbBlock.R` | Yes (S3) | mirrors `print.blmb()` | Combined print for `glmbBlock()` result. |
| `summary.blmb()` | `summary.blmb.R` | Yes (S3) | mirrors `summary.bglmb()` | Combined summary for `lmbBlock()` result. |
| `summary.bglmb()` | `summary.bglmb.R` | Yes (S3) | mirrors `summary.blmb()` | Combined summary for `glmbBlock()` result. |
| `print.summary.blmb()` | `summary.blmb.R` | Yes (S3) | mirrors `print.summary.bglmb()` | Renders `summary.blmb()` output. |
| `print.summary.bglmb()` | `summary.bglmb.R` | Yes (S3) | mirrors `print.summary.blmb()` | Renders `summary.bglmb()` output. |
| `block_check_identifiability()` | `lmbBlock.R` | No | shared by both families | Optional two-level rank preflight, formula interface. |
| `block_check_identifiability_xy()` | `lmbBlock.R` | No | shared by both families | Same, matrix interface. |
| `.blmb_formula_block_meta()` | `lmbBlock.R` | No | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()`, `.blmb_blocks_full_rank()` | Block/formula reconciliation. |
| `.blmb_resolve_block()` | `lmbBlock.R` | No | (called only via the above) | Accepts the 5 `block` shapes. |
| `.blmb_rows_to_data_subset()` | `lmbBlock.R` | No | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()` | Model-frame rows → `data` rows. |
| `.blmb_lmb_display_call()` | `lmbBlock.R` | No | mirrors `.blmb_glmb_display_call()` | Per-block `$call` for `lmb()`. |
| `.blmb_glmb_display_call()` | `glmbBlock.R` | No | mirrors `.blmb_lmb_display_call()` | Per-block `$call` for `glmb()`. |
| `.blmb_assemble()` | `lmbBlock.R` | No | mirrors `.bglmb_assemble()` | Package `k` `lmb()` fits into `"blmb"`. |
| `.bglmb_assemble()` | `glmbBlock.R` | No | mirrors `.blmb_assemble()` | Package `k` `glmb()` fits into `"bglmb"`. |
| `.blmb_coef_means_matrix()` | `lmbBlock.R` | No | `print`/`summary` for both `blmb` and `bglmb` | Stack per-block coefficients. |
| `.blmb_dic_table()` | `lmbBlock.R` | No | `print`/`summary` for both `blmb` and `bglmb` | Stack per-block `DIC`/`pD`. |
| `.blmb_blocks_full_rank()` | `lmbBlock.R` | No | `block_check_identifiability()` | Level-1 rank table, formula interface. |
| `.blmb_blocks_full_rank_xy()` | `lmbBlock.R` | No | `block_check_identifiability_xy()` | Level-1 rank table, matrix interface. |
| `.blmb_resolve_block_calibration_arg()` | `prior_setup_lmebayes.R` | No | `Prior_SetupBlock()` only | Shared-or-per-block resolution for `Prior_Setup()`'s 6 calibration args. |
| `.mrglmb_normalize_pfamily_lists` | `block_core_pfamily.R` (alias; defined in `glmbayesCore`) | No | `lmbBlock()`, `glmbBlock()` | Validate/recycle `pfamily`/`pfamily_list`. |
| `.validate_pfamily_for_rlmb` | `block_core_pfamily.R` (alias; defined in `glmbayesCore`) | No | `lmbBlock()`, `glmbBlock()` | Per-block pfamily validator callback. |
| `lmebayesCore::normalize_block()` | *(package: `lmebayesCore`)* | Yes, there | `lmbBlock()`, `glmbBlock()`, `Prior_SetupBlock()` (via `.blmb_formula_block_meta()`); `block_check_identifiability_xy()` directly | Canonicalize the block partition. Pure R, no C++ — full trace in `lmebayesCore/inst/LMBBLOCK_LMEBAYESCORE_DEPENDENCIES.md`. |

---

## Related files

| Topic | Path |
|-------|------|
| Package-level function tables this expands on | `README.md` ("Row-block models" / "Model setup and mixed-model priors") |
| Exported/internal function inventories | `inst/R_EXPORTED_AND_DOCUMENTED.md`, `inst/R_INTERNAL_HELPERS.md` |
| `lmbBlock()`/`glmbBlock()` source | `R/lmbBlock.R`, `R/glmbBlock.R` |
| `Prior_SetupBlock()` source | `R/prior_setup_lmebayes.R` |
| Shared `pfamily_list` validation source | `R/block_core_pfamily.R` |
| Single `lmebayesCore` dependency, fully traced | `lmebayesCore/inst/LMBBLOCK_LMEBAYESCORE_DEPENDENCIES.md` |
| Two-level identifiability derivation | `lmebayesCore/inst/BLOCK_GIBBS_ERGODICITY.md` |
| `lmb()`/`glmb()`/`Prior_Setup()` (the actual engines) | `glmbayes/R/lmb.R`, `glmbayes/R/glmb.R`, `glmbayesCore/R/Prior_Setup.R` |

---

## Changelog

| Date | Note |
|------|------|
| 2026-07-20 | Revised from a `lmbBlock()`-centric writeup into a symmetric overview covering `lmbBlock()`, `glmbBlock()`, and `Prior_SetupBlock()` together; added the shared-vs-mirrored function classification and `Prior_SetupBlock()`'s own `.blmb_resolve_block_calibration_arg()` helper (previously undocumented in `R_INTERNAL_HELPERS.md`). |
| 2026-07-20 | Initial draft (`lmbBlock()`-centric; superseded by the above). |

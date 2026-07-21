# Manual validation tests

Occasional regression scripts for **major algorithm / sampler changes**. These are **not** part of `R CMD check` or `testthat`.

## Sampler settings

- **No `set.seed()`** — stochastic checks reflect real sampling variability.
- **`lmerb()` / `glmerb()` defaults** for `tv_tol` (0.01), `gap_tol` (0.0196), `mode_gap_max` (1.0).
- Explicit args: **`n`**, **`progbar = TRUE`** only.

| Script / section | Main draws `n` | Route |
|------------------|----------------|-------|
| lmerb §1 Gaussian | 1000 | `rLMMNormal_reg_known_vcov` |
| lmerb §2 ING (same BWC model as §1) | 3000 | `rLMMNormal_reg_estimated_vcov` |
| lmerb dGamma §1 (same BWC model) | 1000 | `rLMMindepNormalGamma_reg_known_vcov` |
| lmerb dGamma §2 (same BWC + ING Block~2) | 3000 | `rLMMindepNormalGamma_reg_estimated_vcov` |
| glmerb §1 Poisson (full airbnb) | 3000 | `rGLMM_reg_known_vcov` |
| glmerb §1 Poisson (`small` → airbnb_small) | 1000 | `rGLMM_reg_known_vcov` |
| glmerb §2 ING (same airbnb) | 3000 / 1000 (`small`) | `rGLMM_reg_estimated_vcov` |

Run from the **lmebayes package root** (reload after code changes):

```r
devtools::load_all("c:/Rpackages/glmbayesCore")   # optional, if testing dev core
devtools::load_all("c:/Rpackages/lmebayes")
setwd("c:/Rpackages/lmebayes")
source("tests/manual/test_lmerb_mer_re_validation.R")
source("tests/manual/test_lmerb_dgamma_mer_re_validation.R")
source("tests/manual/test_glmerb_mer_re_validation.R")
```

Or from PowerShell:

```powershell
cd c:\Rpackages\lmebayes
$env:GLMBAYESCORE_ROOT = "c:\Rpackages\glmbayesCore"

Rscript tests/manual/test_lmerb_mer_re_validation.R
Rscript tests/manual/test_lmerb_dgamma_mer_re_validation.R
Rscript tests/manual/test_glmerb_mer_re_validation.R small   # both sections n=1000
Rscript tests/manual/test_glmerb_mer_re_validation.R           # Poisson n=3000 + ING n=3000
Rscript tests/manual/run_all.R small
```

## Scripts

| Script | glmbayesCore routes | Notes |
|--------|---------------------|--------|
| `test_lmerb_mer_re_validation.R` | known + ING LMM | Same Big Word Club model; §2 ING Block~2 only |
| `test_lmerb_dgamma_mer_re_validation.R` | dGamma LMM + dGamma+ING | Same BWC model; dGamma σ² route only |
| `test_glmerb_mer_re_validation.R` | known + ING GLMM | Block~2 vs ICM + RE cor (Poisson airbnb) |
| `_bwc_lmerb_fixture.R` | (helper) | Shared BWC dat/form; drops rank-deficient `school_id` levels before prior setup |
| `_lmerb_dgamma_fixture.R` | (helper) | dGamma `dispersion_ranef` builder |
| `_block2_fixef_validate.R` | (helper) | Bind `lmebayes:::.validate_manual_block2_fixef` after `load_all()` |
| `_lmerb_re_validate.R` | (helper) | Bind `lmebayes:::.validate_lmerb_re` after `load_all()` |
| `_glmerb_re_validate.R` | (helper) | Bind `lmebayes:::.validate_glmerb_re` after `load_all()` |
| `test_lmerb_fixed_vector_dispersion.R` | (standalone) | Fixed per-group `dispersion_ranef` vector (`dispersion_mode = "fixed_vector"`); `lme4::sleepstudy`, no `bayesrules` needed. Covers `lmerb()`, `summary_sigma2()`, validation errors, and `glmerb(family = gaussian())` parity |

Demos under `demo/` mirror these routes (`Ex_12`/`Ex_21` Gaussian+ING; `Ex_24`/`Ex_25` dGamma) with fuller printed output.

**BWC fixture:** `_bwc_lmerb_fixture.R` subsets to algebraically **full-rank** schools (`model_setup()$re_rank`) before `Prior_Setup_lmebayes`, matching `Ex_24` / `Ex_25`.

**dGamma note:** Block~1 uses center → build → dispersion build → sim; the manual dGamma script uses the **same** filtered Big Word Club rows and formula as `test_lmerb_mer_re_validation.R`.

**dGamma §1–§2 (BlockEnvelopeCentering):** chain means should match **`lmer_full`** (cor ≈ 1). **`fixef.mode` / `ranef.mode` ICM checks are skipped** (`z_icm_max = Inf`, relaxed RE ICM). **`se_ratio` uses the full `[0.85, 1.05]` band on all mapped rows** (do not pass `ing = TRUE` here — that only relaxes `null_effects` for Gaussian/glmer ING shrinkage). §2 still differs from §1 on dispersion sampling and structural checks only.

**Block~2 fixef `se_ratio`:** Gaussian §1 and glmer §1 enforce `0.85 <= post_sd/lmer_se <= 1.05` on all mapped rows. **Gaussian/glmer ING** (`ing = TRUE`) relaxes the lower bound on `null_effects` only. dGamma §1–§2 skip ICM via `z_icm_max = Inf` but keep full `se_ratio` enforcement.

## Environment

| Variable | Purpose |
|----------|---------|
| `LMEBAYES_ROOT` | Package root if not running from it |
| `GLMBAYESCORE_ROOT` | Load dev glmbayesCore before lmebayes (see `_load.R`) |

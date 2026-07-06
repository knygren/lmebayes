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
| `_bwc_lmerb_fixture.R` | (helper) | Shared Big Word Club dat/form for lmerb manual tests |
| `_lmerb_dgamma_fixture.R` | (helper) | dGamma `dispersion_ranef` builder |
| `_glmerb_re_validate.R` | (helper) | glmerb ordering + `glmer_full` cor |

Demos under `demo/` mirror these routes (`Ex_12`/`Ex_21` Gaussian+ING; `Ex_24`/`Ex_25` dGamma) with fuller printed output.

**dGamma note:** Block~1 still uses **BlockEnvelopeCentering** until BlockEnvelopeSim ships; the manual dGamma script now uses the **same** Big Word Club rows and formula as `test_lmerb_mer_re_validation.R` so results are directly comparable across routes.

## Environment

| Variable | Purpose |
|----------|---------|
| `LMEBAYES_ROOT` | Package root if not running from it |
| `GLMBAYESCORE_ROOT` | Load dev glmbayesCore before lmebayes (see `_load.R`) |

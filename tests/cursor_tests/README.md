t# cursor_tests — manual regression via package demos

This folder is **not** part of `R CMD check`. Do **not** add `test-*.R`, smoke
scripts, or ad-hoc fixtures here.

## Policy

All manual regression for **`lmerb()`**, **`glmerb()`**, and the six underlying
**glmbayesCore** reg-route engines must use the demos and help examples listed
in [README.md — Examples and Demos](../README.md#examples-and-demos), plus
**`Ex_23`** for the full Gaussian LMM 2×2. Standard workflow:

```r
ps <- Prior_Setup_GLMM(form, data = dat, ...)
pf <- pfamily_list(ps)   # or ptypes = "dIndependent_Normal_Gamma" for ING Block~2
fit <- lmerb(form, data = dat, pfamily_list = pf, dispersion_ranef = ps$dispersion_ranef, ...)
# or glmerb(..., family = ..., pfamily_list = pf, ...)
```

## Six **glmbayesCore** routes ↔ demos

| Matrix route (in **glmbayesCore**) | What you set in **lmebayes** | Primary demo(s) |
|------------------------------------|------------------------------|-----------------|
| `rLMMNormal_reg_known_vcov()` | Gaussian; `pfamily_list(ps)`; scalar `dispersion_ranef` | `example("lmerb")`; `Ex_14_lmerb_Sleepstudy`; `Ex_12_lmerb_BigWordClub`; **Ex_23 case 1** |
| `rLMMNormal_reg_estimated_vcov()` | Gaussian; ING Block~2 (`pfamily_list(..., ptypes = "dIndependent_Normal_Gamma")`); scalar σ² | `Ex_20_lmerb_ING_pilot`; `Ex_21_lmerb_ING_BigWordClub`; **Ex_23 case 2** |
| `rLMMindepNormalGamma_reg_known_vcov()` | Gaussian; dGamma `dispersion_ranef`; dNormal Block~2 | **Ex_23 case 3** (`simulate = FALSE`); Gibbs: same call with `simulate = TRUE` |
| `rLMMindepNormalGamma_reg_estimated_vcov()` | Gaussian; dGamma σ² + ING Block~2 | **Ex_23 case 4** (`simulate = FALSE`); Gibbs: same with `simulate = TRUE` |
| `rGLMM_reg_known_vcov()` | Non-Gaussian; `pfamily_list(ps)` (dNormal Block~2) | `example("glmerb")`; `Ex_14_glmerb_airbnb_small`; `Ex_16_glmerb_book_banning` |
| `rGLMM_reg_estimated_vcov()` | Non-Gaussian; ING Block~2 | `Ex_22_glmerb_book_banning_ING` |

**Ex_23** = `demo("Ex_23_lmerb_joint_posterior_mode_four_cases", package = "lmebayes")`
covers all four LMM σ² × Block~2 combinations at the **`lmerb()`** API (joint
posterior mode / ICM). It is the canonical reference for the two **dGamma σ²**
routes; add `simulate = TRUE` when validating full Gibbs sampling.

## How to run

```r
library(lmebayes)

## Quick (help examples; no long Gibbs)
example("lmerb")
example("glmerb")

## Full workflows (see README for the full list)
demo("Ex_14_lmerb_Sleepstudy", package = "lmebayes")
demo("Ex_16_glmerb_book_banning", package = "lmebayes")
demo("Ex_20_lmerb_ING_pilot", package = "lmebayes")
demo("Ex_22_glmerb_book_banning_ING", package = "lmebayes")
demo("Ex_23_lmerb_joint_posterior_mode_four_cases", package = "lmebayes")
```

After sampling runs, inspect convergence with `print(fit$sweep_history$main)` and
`plot_mean_convergence()` / `plot_var_convergence()` (see README demos **Ex_16**, **Ex_21**, **Ex_22**).

## Official **testthat** suite

Automated checks live under **`tests/testthat/`** only (fast unit tests:
`test-sweep-history.R`, `test-summary-lmerb.R`, `test-prior-setup-*.R`, etc.).
Long Gibbs / MER validation is under **`tests/manual/`** (not run by `devtools::check()`).
See [tests/manual/README.md](../manual/README.md). Quick start:

```r
## From package root (PowerShell)
Rscript tests/manual/test_lmerb_mer_re_validation.R
Rscript tests/manual/test_lmerb_dgamma_mer_re_validation.R
Rscript tests/manual/test_glmerb_mer_re_validation.R small
Rscript tests/manual/run_all.R small
```

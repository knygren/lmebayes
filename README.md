# lmebayes

![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/lmebayes?label=version)
![R-universe](https://knygren.r-universe.dev/badges/lmebayes)
![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)

**lmebayes** provides near-independent posterior samples for Bayesian linear and generalized linear
**mixed-effects** models via **two-block Gibbs sampling** (sampling engines in **glmbayesCore**).
Its primary interfaces, `lmerb()` and `glmerb()`, are Bayesian analogues of **lme4** `lmer()` and
`glmer()`, supporting Gaussian, Poisson, binomial, and Gamma response families under log-concave
likelihoods. Row-block BY-style fits use `lmbBlock()` and `glmbBlock()`; matrix-level block samplers
include `rNormalRegBlock()` and `rNormalGLMBlock()`.

Priors, `pfamily` objects, and iid GLM sampling within blocks come from the **glmbayes** dependency.
Mixed-model methodology and background vignettes are in **glmbayes** (Chapters 17 and 18 for LMMs and
GLMMs). **lmebayes** does not ship vignettes yet; use function help, this README, and the package
demos. For Gaussian models, inner Gibbs sweep counts can be calibrated from a total-variation
tolerance (`tv_tol`); non-Gaussian GLMMs may run a pilot stage when `gap_tol` is set (see `?glmerb`).
Non-Gaussian **`glmerb()`** / **`rglmerb()`** use **`glmbayesCore::rGLMM()`** (sweep-outer
driver) only; legacy **`rglmerb_v5()`** / C++ short-chain code has been removed (see [NEWS.md](NEWS.md)).

This repository is **0.1.0** in development.
The [GitHub](https://github.com/knygren/lmebayes) repository holds the source;
[R-Universe](https://knygren.r-universe.dev/lmebayes) builds binaries from it.
See [NEWS.md](https://github.com/knygren/lmebayes/blob/main/NEWS.md) for changes.

## Function overview

The tables below list symbols exported from **lmebayes** (see `help(package = "lmebayes")`).
They follow the same broad grouping used in **glmbayes** vignette
[Chapter A01](https://knygren.r-universe.dev/articles/glmbayes/Chapter-A01.html):
core formula interfaces, prior helpers, low-level samplers, diagnostics, and
advanced simulation callbacks.

### Functions defined in **lmebayes**

#### Mixed-effects model fitting (`lme4`-style)

| Function | Role |
|----------|------|
| `lmerb()` | Bayesian linear mixed-effects model (LMM); two-block Gibbs sampler. Bayesian analogue of `lme4::lmer()`. |
| `glmerb()` | Bayesian generalized linear mixed-effects model (GLMM). Bayesian analogue of `lme4::glmer()`. |
| `print()` / `summary()` | S3 methods for `"lmerb"` and `"glmerb"` objects (posterior summaries, sweep history). |

Set `simulate = FALSE` on `lmerb()` / `glmerb()` for ICM posterior mean/mode only (no stored draws).

#### Row-block models (SAS `BY`-style splits)

| Function | Role |
|----------|------|
| `lmbBlock()` | One `lmb()` fit per row block (shared formula, block-specific priors). |
| `glmbBlock()` | One `glmb()` fit per row block. |
| `print()` / `summary()` | S3 methods for `"blmb"` and `"bglmb"` (lists of block fits). |

#### Model setup and mixed-model priors

| Function | Role |
|----------|------|
| `model_setup()` | Parse an `lme4`-style formula into design matrices, variance components, and identifiability checks (single grouping factor). |
| `Prior_Setup_lmebayes()` | Calibrate Block~2 hyperpriors and observation-level dispersion from a reference `lmer` / `glmer` fit. |
| `Prior_SetupBlock()` | Run `Prior_Setup()` independently on each row block (for `lmbBlock()` / `glmbBlock()`). |
| `pfamily_list()` | S3 method for `"lmebayes_prior_setup"`: build per–random-effect `pfamily` objects from `Prior_Setup_lmebayes()`. |
| `print()` | S3 methods for `"model_setup"` and `"lmebayes_prior_setup"`. |

Typical workflow: `model_setup()` → `Prior_Setup_lmebayes()` → `pfamily_list(ps)` → `lmerb()` / `glmerb()`.

#### Low-level mixed-effects samplers (matrix / design API)

These mirror `rlmb()` / `rglmb()` in **glmbayes**: minimal overhead for repeated Gibbs calls.
`lmerb()` and `glmerb()` invoke them after formula parsing and prior construction.

| Function | Role |
|----------|------|
| `rlmerb()` | Gaussian LMM sampler (two-block Gibbs; replicate chains). |
| `rglmerb()` | GLMM sampler: Gaussian → `rLMMNormal_reg()` / `rLMMindepNormalGamma_reg()`; other families → `glmbayesCore::rGLMM()` (sweep-outer; optional pilot when `gap_tol` is set). |

#### Diagnostics and build utilities

| Function | Role |
|----------|------|
| `plot_sweep_history_diag()` | Histogram diagnostics on inner Gibbs sweep history from `lmerb()` / `glmerb()`. |
| `has_opencl()` | Whether **this** **lmebayes** build was compiled with OpenCL (distinct from runtime GPU probes in **opencltools**). |

### Re-exported functions

Symbols below are implemented in **glmbayes** or **glmbayesCore** and re-exported so a single
`library(lmebayes)` load covers mixed models and the iid GLM tools they build on.

#### From **glmbayes** — iid Bayesian `lm` / `glm`

| Function | Role |
|----------|------|
| `lmb()` | Bayesian linear model (iid draws). Analogue of `lm()`. |
| `glmb()` | Bayesian GLM (iid draws). Analogue of `glm()`. |
| `directional_tail()` | Directional tail diagnostic for prior–posterior disagreement (see **glmbayes** Chapter A04). |

Row-block wrappers `lmbBlock()` and `glmbBlock()` call these per block.

#### From **glmbayesCore** — prior families and default calibration

| Function | Role |
|----------|------|
| `Prior_Setup()` | Default prior calibration for a GLM/LM formula (Zellner-style `mu`, `Sigma`, dispersion, conjugate components). |
| `dNormal()`, `dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, `dGamma()` | `pfamily` constructors passed to `lmb()`, `glmb()`, and block samplers. |
| `pfamily_list()` | Generic (from **glmbayesCore**); **lmebayes** adds the `lmebayes_prior_setup` method above. |

See **glmbayes** README sections *Supported families, links, and pfamilies* and *Prior_Setup* for wiring details.

For the full simulation and envelope map, see **glmbayes** vignettes
[Chapter A05](https://knygren.r-universe.dev/articles/glmbayes/Chapter-A05.html) and
[Chapter A08](https://knygren.r-universe.dev/articles/glmbayes/Chapter-A08.html).

## Installation

**CRAN (release 0.9.5)**

```r
install.packages("glmbayes")
```

**GitHub / R-Universe** (install from both CRAN and R-Universe repositories if you want R-Universe binaries or faster mirrors):

```r
install.packages("glmbayes",
                 repos = c("https://cloud.r-project.org",
                           "https://knygren.r-universe.dev"))
```

Prebuilt binaries from CRAN (0.9.5) and R-Universe are built **without OpenCL GPU
support**. For the CRAN release, OpenCL requires installing **from source** on a
system with OpenCL development files available. To set up GPU acceleration, follow

**Chapter 16 — Large models: GPU acceleration using OpenCL**
https://knygren.r-universe.dev/articles/glmbayes/Chapter-16.html

## Minimal Working Example

Requires the **bayesrules** package (`install.packages("bayesrules")`).

    library(lmebayes)

    data(big_word_club, package = "bayesrules")
    dat <- subset(
      big_word_club,
      !is.na(score_ppvt) & !is.na(invalid_ppvt) & invalid_ppvt == 0L
    )
    dat$school_id <- factor(dat$school_id)
    dat <- dat[complete.cases(dat[, c("score_ppvt", "distracted_ppvt",
                                      "free_reduced_lunch", "school_id")]), ]

    form <- score_ppvt ~ free_reduced_lunch + distracted_ppvt +
      (1 + distracted_ppvt || school_id)

    ## Classical lmer (reference fit embedded in lmerb)
    lme4::lmer(form, data = dat)

    ## Bayesian lmerb — prior setup + ICM posterior mean/mode (no Gibbs draws)
    ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
    fit <- lmerb(
      form,
      data             = dat,
      pfamily_list     = pfamily_list(ps),
      dispersion_ranef = ps$dispersion_ranef,
      simulate         = FALSE
    )

    lmebayes:::print_coef_means(fit)
    print(fit)
    summary(fit)

`Prior_Setup_lmebayes()` calibrates Block~2 hyperpriors from a weak-prior **lmer** fit;
`lmerb(..., simulate = FALSE)` returns that reference fit plus exact **ICM** posterior
mean/mode values (no stored draws). For iid Gibbs samples, set `simulate = TRUE` or run
the demos listed below.

## Priors and GLM families (`glmbayes`)

Formula-based priors (`Prior_Setup`, `pfamily`, `dNormal`, etc.) and `glmb()` / `lmb()` live in the
**glmbayes** dependency. **lmebayes** adds row-block priors via `block_prior_setup()` and block Gibbs
samplers. See `?glmbayes::Prior_Setup`, `?glmbayes::pfamily`, and `vignette("Chapter-04", package = "glmbayes")`.


## Examples and Demos

Use `example()` for quick help-page examples (ICM / setup only; safe for `R CMD check`).
Use `demo()` for full Gibbs workflows with stored draws (may take minutes).

    ## Bayesian Linear Mixed-effects model (no Simulation) 

    example("lmerb")    ## big_word_club Gaussian LMM (small formula)

    ## Small lme4-style model with simulation (sleepstudy)

    demo("Ex_14_lmerb_Sleepstudy", package = "lmebayes")

    ## Bayesian Generalized Linear Mixed-effects model (no Simulation) 

    example("glmerb")   ## airbnb_small Poisson GLMM

    ## Same model with simulation
  
    demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")

    ## Larger lmerb model

    demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")

    ## Larger glmerb model

    demo("Ex_13_glmerb_Airbnb", package = "lmebayes")

    ## Bayes Rules book_banning (binomial GLMM; requires bayesrules)

    demo("Ex_16_glmerb_book_banning", package = "lmebayes")
    ## Minimal Ch. 18-style: violent_i + state RE / random slope

    demo("Ex_19_glmerb_book_banning_state_covariates", package = "lmebayes")
    ## State covariates in X_hyper for the intercept RE

    ## Independent Normal–Gamma Block~2 priors (ING; pilot + main sampling)

    demo("Ex_20_lmerb_ING_pilot", package = "lmebayes")
    ## Gaussian LMM, ING + pilot (small Big Word Club model)

    demo("Ex_21_lmerb_ING_BigWordClub", package = "lmebayes")
    ## Gaussian LMM, ING + pilot (full Big Word Club; compare Ex_12)

    demo("Ex_22_glmerb_book_banning_ING", package = "lmebayes")
    ## Binomial GLMM, ING on RE components (compare Ex_16)


## Methodology

For generalized linear models where well known sampling methods are unavailable, sampling follows the
framework from Nygren and Nygren (2006), using likelihood subgradients to construct enveloping functions for
the posterior distribution. When the posterior is approximately normal, the expected number of draws per
acceptance is bounded as per that paper and as discussed in the **glmbayes** vignettes.
Dispersion can be sampled via `rGamma_reg()` (standalone) or jointly with coefficients via
`rNormalGamma_reg()` and `rindepNormalGamma_reg()`.

## GPU Acceleration Using OpenCL

The implemented algorithms tend to have acceptable performance on CPUs up to around 10-14 dimensions.
For larger models, the envelope construction is embarrassingly parallel. To accelerate envelope construction
in such cases, the package provides optional GPU acceleration using OpenCL. This requires that users have
GPU enabled machines and an OpenCL installation. See `vignette("Chapter-16", package = "glmbayes")`
and `vignette("Chapter-A10", package = "glmbayes")` in **glmbayes**.

## Documentation

**lmebayes** does not ship vignettes yet; use the [Function overview](#function-overview) above,
`?lmebayes`, and function help pages here.
For GLM/Gibbs sampler background and tutorials, see **glmbayes**: `browseVignettes("glmbayes")`
or https://knygren.r-universe.dev/articles/glmbayes/index.html .
Mixed-model methodology is covered in **glmbayes** Chapters 17 (LMMs) and 18 (GLMMs).

**Maintainers:** `R/` symbols are split into
[inst/R_EXPORTED_AND_DOCUMENTED.md](inst/R_EXPORTED_AND_DOCUMENTED.md) (exports and
`man/` pages) and
[inst/R_INTERNAL_HELPERS.md](inst/R_INTERNAL_HELPERS.md) (`@noRd` / undocumented
helpers). Index: [inst/R_FUNCTION_INVENTORY.md](inst/R_FUNCTION_INVENTORY.md).

## Feature Highlights

- S3 interface mirroring the structure of base glm()
- Posterior predictive checks via `pp_check()` from the 'bayesplot' package for fitted `glmb` objects
- Accept-reject sampling for log-concave likelihoods
- Samplers for both fixed and variable dispersion
- Reuses **glmbayes** samplers and vignetted methodology (mixed-model vignettes planned for **lmebayes**)
- Modular prior setup function

## Limitations

- Non-log-concave likelihoods are not currently supported

## Future Plans

- **R Mathlib (`nmath`) usage from C:** Today the package vendors local copies of
  selected R Mathlib routines and headers in `*.c` sources. The plan is to switch
  to calling the **same `nmath` functions that ship with R**, via the supported
  linking/API path, so maintenance tracks base R instead of duplicating sources.
- **OpenCL / GPU code upstream:** Routines currently living under the
  **openclport** and **nmathopencl** namespaces are slated to move into dedicated
  upstream packages. **nmathopencl** is already available on
  [R-Universe](https://knygren.r-universe.dev/nmathopencl); a **CRAN** release is targeted,
  after which glmbayes can depend on that package for a substantial share of
  OpenCL- and GPU-related functionality rather than carrying those implementations
  here.
- **Conjugate priors for intercept-only GLMs:** Add **pfamily** specifications
  that supply conjugate priors for **intercept-only** `glm()`-style models (a
  single mean structure / scalar linear predictor), complementing the existing
  prior families for general designs.
- **bayestestR integration:** Add methods or small wrappers so **bayestestR**
  summaries and diagnostics can be used with **`glmb` / `lmb`** fits in the same
  way as with other Bayesian modeling workflows.

Further performance and algorithm work:

- Poisson speed (OpenCL and simulation): Precompute the log-factorial term `log(y!)`
  once per observation and reuse it in both OpenCL envelope construction and
  accept-reject simulation, since it depends only on the response, to reduce
  redundant `lgamma` evaluation and improve performance for large Poisson models.
- Grid selection (simulation): Precompute cumulative PLSD and use inverse CDF
  sampling (e.g. binary search) to select the grid component per candidate
  instead of scanning PLSD, improving the simulation loop when many candidates
  are evaluated.
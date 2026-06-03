# glmbayes

![CRAN status](https://www.r-pkg.org/badges/version/glmbayes)
![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/glmbayes)
![Monthly downloads](https://cranlogs.r-pkg.org/badges/glmbayes)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayes?label=version)
![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayes/R-CMD-check.yaml?label=R%20CMD%20Check)

glmbayes provides independent and identically distributed (iid) samples for Bayesian Generalized Linear Models (GLMs).
Its primary interface, glmb(), serves as a Bayesian analogue to R's glm() function, supporting Gaussian, Poisson,
Binomial, and Gamma families under log-concave likelihoods. Sampling for most models is performed using accept-reject
methods based on likelihood subgradients (Nygren and Nygren, 2006). For Gaussian models, the package also includes
lmb(), a Bayesian counterpart to R's lm().

The package includes a rich set of supporting tools for prior specification, model diagnostics, and method functions
that mirror those for lm() and glm(). Most functions are extensively documented.
Background vignettes for the underlying samplers live in the **glmbayes** package (`browseVignettes("glmbayes")`).
**lmebayes** vignettes are planned separately.

This repository is **0.9.6** in development. The current **CRAN release is version 0.9.5**
([CRAN](https://CRAN.R-project.org/package=glmbayes)).
The [GitHub](https://github.com/knygren/glmbayes) repository holds the source; [R-Universe](https://knygren.r-universe.dev/glmbayes) builds binaries from it.
See [NEWS.md](https://github.com/knygren/glmbayes/blob/main/NEWS.md) for changes.

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

    library(glmbayes)

    # Dobson (1990), p. 93: Randomized Controlled Trial
    counts <- c(18,17,15,20,10,20,25,13,12)
    outcome <- gl(3,1,9)
    treatment <- gl(3,3)
    print(d.AD <- data.frame(treatment, outcome, counts))

    ## Classical glm
    glm.D93 <- glm(counts ~ outcome + treatment,
                   family = poisson())

    ## Bayesian glmb (via glmbayes)
    ps <- glmbayes::Prior_Setup(counts ~ outcome + treatment, family = poisson())
    glmb.D93 <- glmbayes::glmb(counts ~ outcome + treatment,
                     family = poisson(),
                     pfamily = glmbayes::dNormal(mu = ps$mu, Sigma = ps$Sigma))

    summary(glmb.D93)

## Priors and GLM families (`glmbayes`)

Formula-based priors (`Prior_Setup`, `pfamily`, `dNormal`, etc.) and `glmb()` / `lmb()` live in the
**glmbayes** dependency. **lmebayes** adds row-block priors via `block_prior_setup()` and block Gibbs
samplers. See `?glmbayes::Prior_Setup`, `?glmbayes::pfamily`, and `vignette("Chapter-04", package = "glmbayes")`.


## Examples and Demos

Use `example()` and `demo()` to explore built-in examples and demos for supported families and links:

    ## Bayesian linear regression
    example("lmb")

    ## Bayesian generalized linear models
    example("glmb")

    ## Predictions for fitted glmb objects (newdata, type, etc.)
    example("predict.glmb")

    ## Deviance residuals and simulate() for posterior predictive checks (menarche)
    example("residuals.glmb")

    ## Two-block Gibbs sampler compared with iid sampling (linear model)
    example("rlmb")

    ## Default prior specification (glmbayes)
    example("Prior_Setup", package = "glmbayes")

    ## Matrix-input GLM example with an informative prior
    example("rglmb")

    ## Two-step Boston example: estimates and summarizes models with unknown
    ## dispersion using dGamma priors via rGamma_reg, rglmb, rlmb, glmb, and lmb
    example("summary.rGamma_reg")

    ## High-dimensional Gaussian model (14 predictors) with GPU acceleration (requires OpenCL)
    example("Boston_centered")

    ## High-dimensional binomial model (14 predictors) with GPU acceleration (requires OpenCL)
    example("Cleveland")

    ## Hierarchical linear model (Rubin/Gelman 8-schools) via rlmb
    demo("Ex_07_Schools")

    ## Hierarchical generalized linear model (Poisson BikeSharing) via rglmb
    demo("Ex_09_BikeSharingPoisson")

    ## Detailed simulation pipeline for rNormalGLM models (JASA 2006; glmbayes Chapter A05)
    example("rNormalGLM_std")

    ## Detailed simulation pipeline for rIndepNormalGammaReg models (glmbayes Chapter A07)
    example("rIndepNormalGammaReg_std")

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

**lmebayes** does not ship vignettes yet; use `?lmebayes` and function help pages here.
For GLM/Gibbs sampler background and tutorials, see **glmbayes**: `browseVignettes("glmbayes")`
or https://knygren.r-universe.dev/articles/glmbayes/index.html .

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
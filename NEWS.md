# lmebayes (development version)

* **Per-component `pwt` and decoupled dispersion prior in
  `Prior_Setup_lmebayes()`:** `pwt` now accepts, besides a scalar, a list
  with one element per random-effect component (named with the RE
  coefficient names or positional); each element is a scalar (recycled
  over that component's Block-2 predictors) or a per-predictor vector
  (optionally named with the `X_hyper[[k]]` columns).  `Sigma_fixef` is
  scaled elementwise by `sqrt((1-w_i)/w_i) * sqrt((1-w_j)/w_j)`, matching
  the `glmbayesCore::Prior_Setup` vector-`pwt` convention and reducing to
  the classic `(1-pwt)/pwt` factor for equal weights.  Two new optional
  arguments decouple the Block-2 dispersion (precision) prior from the
  coefficient weights: `pwt_dispersion` (relative weight in (0,1)) and
  `n_prior_dispersion` (absolute effective prior sample size in group
  units), each a scalar or per-component list/vector; at most one may be
  supplied, and when neither is, the values are derived from the
  per-component mean `pwt` for consistency.  The returned object always
  carries mutually consistent `$pwt_dispersion` and `$n_prior_dispersion`
  per-component vectors (`n_k = J w_k / (1 - w_k)`), which
  `pfamily_list()` now uses to calibrate `dIndependent_Normal_Gamma`
  shape/rate per component instead of re-deriving them from the scalar
  `pwt`.  `print()` shows per-component weights and the dispersion-prior
  source.  Covered by `data-raw/test_prior_setup_pwt.R`.

* **Conservative TV calibration for `dIndependent_Normal_Gamma` priors:**
  `lmerb()`/`glmerb()` now accept ING components in `pfamily_list`,
  provided each supplies a positive `disp_lower` (lower truncation of the
  dispersion).  `disp_lower` replaces the `dNormal` dispersion as the
  plug-in `tau^2_k` in the eigenvalue / TV calibration: smaller `tau^2`
  increases the block coupling and hence `lambda*`, so the disp_lower-based
  rate upper-bounds the contraction rate for every dispersion in the
  truncated support.  Because Block-2 dispersion sampling is not
  implemented yet, the fit displays the calibration
  (`conservative: ING tau^2_k = disp_lower`) and stops, returning the ICM
  mode plus `$convergence` (method `"disp_lower_bound"`, or
  `"<base>+disp_lower_bound"` in `glmerb`) without draws.  On
  `big_word_club` with `disp_lower = tau^2_k / 2`, `lambda*` rises from
  0.839 to 0.903 and `m_min` from 11 to 18, matching an explicit `dNormal`
  fit at `tau^2/2` exactly (`data-raw/test_ing_calibration.R`).
  `pfamily_list()` now fills in a default
  `disp_lower = 1 / qgamma(0.99, shape, rate)` - the 0.01 quantile of the
  implied inverse-Gamma dispersion prior (reciprocal of the 99th percentile
  of the Gamma precision prior) - so ING lists it builds pass
  `lmerb()`/`glmerb()` validation out of the box and the calibration covers
  99% of the prior dispersion mass.  Under the diffuse default calibration
  (`pwt = 0.01`) this quantile sits at roughly 3-5% of `tau^2_k` on
  `big_word_club`, giving a strongly conservative `lambda* = 0.989`,
  `m_min = 156`.

* **`lmerb()`/`glmerb()` prior interface migrated to pfamily lists:** The
  `measurement_prior_list` argument (a whole `Prior_Setup_lmebayes()`
  object) is replaced by two explicit arguments placed before `n`:
  `pfamily_list` (named list of `dNormal` pfamilies, one per random-effect
  coefficient - the Block-2 hyperpriors) and `dispersion_ranef` (the
  observation-level measurement dispersion, a known constant for now;
  required for `gaussian()`, must be `NULL` for `poisson()`/`binomial()`).
  The Block-1 random-effect covariance is reconstructed from the pfamily
  dispersions (`Sigma_ranef = diag(tau^2_k)`), so the pair is
  information-complete.  Typical workflow:
  `ps <- Prior_Setup_lmebayes(...)`, then
  `lmerb(f, dat, pfamily_list = pfamily_list(ps),
  dispersion_ranef = ps$dispersion_ranef)`.
  `dIndependent_Normal_Gamma` components are rejected with a clear message
  until Block-2 dispersion sampling is implemented.  The fitted object's
  `$prior` now stores the normalized container (`pfamily_list`,
  `dispersion_ranef`, `Sigma_ranef`, `prior_list`); `summary()` methods are
  unchanged.  Validation/conversion lives in the internal helper
  `.lmebayes_priors_from_pfamily_list()`.

* **Block-2 hyperpriors as pfamily objects (`pfamily_list()`):** New S3
  method `pfamily_list.lmebayes_prior_setup()` converts the per-component
  Block-2 hyperprior parameters of a `Prior_Setup_lmebayes()` object into a
  named list of `glmbayesCore` pfamily objects (one per random-effect
  coefficient).  The `ptypes` argument is either a single string recycled
  to every component or a character vector / list with one entry per
  component (optionally named, in any order); allowed values are
  `"dNormal"` (known Block-2 dispersion `tau^2_k`) and
  `"dIndependent_Normal_Gamma"` (Gamma prior on the Block-2 precision,
  calibrated with the `shape_ING` convention from
  `glmbayesCore::Prior_Setup()` using `n0 = J * pwt/(1-pwt)`:
  `shape = (n0+1)/2 + p_k/2`, `rate = tau^2_k * n0/2`).  The generic lives
  in `glmbayesCore` and is re-exported.

* **Approximate TV calibration for non-Gaussian `glmerb` + `m_convergence`
  override:** Non-Gaussian `glmerb()` now derives its sweep count from the
  same Theorem 3 machinery applied to the *local-Gaussian approximation of
  the posterior at its mode*: per-observation likelihood precisions are
  evaluated at the ICM posterior mode
  (`glmbayesCore::two_block_mode_weights()`) and fed to
  `glmbayesCore::two_block_rate(weights = )`.  The derived `m_min` is the
  minimum number of inner Gibbs sweeps required to converge to that
  hypothetical multivariate normal approximation - a lower bound for the
  true (non-normal) posterior, replacing the previous fixed `10L`.  Both
  `lmerb()` and `glmerb()` gain an optional `m_convergence` argument: when
  supplied it overrides the derived value but is floored at `m_min`
  (`max(m_convergence, m_min)`, warning if raised) - typical use is picking
  a larger number, e.g. double the derived bound.  The calibration is
  reported in a clearly labeled line (`exact` vs
  `approximate (local-Gaussian at mode, <family>)`) and stored in the fitted
  object as `$convergence` (`method`, `tv_tol`, `lambda_star`,
  `eigenvalues`, `m_min`, `m_convergence`).  On `bayesrules::airbnb`
  (Poisson) the heuristic gives `lambda* = 0.48`, `m_min = 4`.

* **TV-calibrated Gibbs sweeps (`tv_tol`):** `lmerb()` and `glmerb()` gain a
  `tv_tol` argument (default `0.01`, the conventional threshold of the
  honest-burn-in literature; Jones & Hobert 2001).  For `lmerb()` (and
  `glmerb()` with `family = gaussian()`) the number of inner Gibbs sweeps per
  stored draw (`m_convergence`, previously hardcoded to `10L`) is now derived
  exactly: the Remark 8 eigenvalue spectrum (Nygren 2020) is computed with
  `glmbayesCore::two_block_rate()` and the exact Theorem 3 TV bound is
  inverted with `glmbayesCore::two_block_l_for_tv()`, plus one sweep for the
  half-step lag of the stored random-effect draw.  Chains start at the exact
  joint posterior mean (ICM), so the bound's mean term vanishes and each
  stored draw is guaranteed within `tv_tol` of the exact joint posterior in
  total variation.  For non-Gaussian `glmerb()` families no exact calibration
  exists; `tv_tol` is accepted but currently ignored (with a message) pending
  an approximate local-Gaussian (IRLS-weight) calibration.  Regression test:
  `data-raw/test_tv_tol_arg.R`.

# glmbayes 0.9.6

## Highlights

* **Row-block (`block_*`) and multi-response (`multi_*`) APIs:** **`block_prior_setup()`**
  and **`block_lmb()`** fit separate **`lmb()`** models per observation block
  (SAS `BY`-style row splits; class **`blmb`**). **`multi_lmb()`** fits several
  response columns with a shared formula (class **`mlmb`**). Gibbs block samplers
  are **`block_rNormalGLM()`** / **`block_rNormalGLM_update()`** (aliases
  **`rNormalGLM_reg_block*`** retained).

* **Conjugate GLM priors (Poisson, binomial, Gamma):** New closed-form IID
  sampling paths for intercept-only models with identity links. **`dBeta()`**
  with **`rBeta_reg()`** supports Beta–Binomial(identity) conjugate updates;
  **`dGamma(Inv_Dispersion = FALSE)`** with **`rGamma_Conjugate_reg()`**
  supports Gamma–Poisson(identity) and Gamma–Gamma(identity) rate priors.
  **`Prior_Setup()`** can calibrate conjugate hyperparameters for these
  families (weighted Poisson rate and binomial probability defaults). See
  **`?dBeta`**, **`?dGamma`**, and the Chapter 02 / Chapter 07–11 vignettes.

* **Vignette structure:** Reworked **Chapter 00** as a roadmap across five
  main parts plus technical appendices. **Chapter 02** is now a conceptual
  introduction to single-parameter conjugacy; worked examples move to
  **Chapter 02-S01** through **Chapter 02-S05** (Beta–Binomial, Normal–Normal,
  Gamma–Poisson, exposure-weighted Poisson, and related topics). A **Companion
  textbooks** section in Chapter 00 indexes optional Bayes Rules! and `LearnBayes`
  appendices tied to the main GLM chapters.

* **`opencltools` import:** Core host/runtime OpenCL discovery and diagnostics
  (`detect_*`, PATH helpers, environment checks) now live in the **`opencltools`**
  package (`Imports`, >= 0.8.0). **glmbayes** keeps package-specific entry
  points (`has_opencl()`, `diagnose_glmbayes()`) that report compile-time
  OpenCL status for this build while delegating shared GPU/runtime checks—reducing
  duplicated maintenance in **glmbayes**.

* **Bayes Rules! companion examples:** Optional vignette appendices reproduce
  book datasets and published posterior summaries using **`lmb()`**, **`glmb()`**,
  **`Prior_Setup()`**, and **`dNormal()`** (suggested package **`bayesrules`** for
  data only). Coverage includes **`bikes`** (Ch. 03), **`weather_perth`** (Ch. 08–09),
  **`equality_index`** (Ch. 10), Gamma–Poisson conjugacy (Ch. 02-S04), and a
  scope note for Gamma regression (Ch. 11). Comparison tables use **printed book
  values**, not live **`rstanarm`** fits. See **Chapter 00** § Companion textbooks.

* **`LearnBayes` examples:** **Chapter 02-S04**, Appendix A, maps the
  **`hearttransplants`** example from Albert (2009) / `LearnBayes` (exposure-weighted
  Gamma–Poisson conjugacy) to **`glmb()`** with analytic Albert posteriors for
  verification (suggested package **`LearnBayes`**).

## Other changes

* Expanded **testthat** coverage for **`dBeta()`** / binomial(identity) conjugate
  paths and related **`glmb()`** integration.

# glmbayes 0.9.5

* **Tests / CRAN:** All **OpenCL**-specific **testthat** blocks now call
  **`skip_on_cran()`** (in addition to **`skip_if_no_opencl()`**), consistent
  with existing Boston/Cleveland OpenCL tests. OpenCL coverage remains for local
  checks and source builds with OpenCL; CRAN checks avoid parallel/GPU-heavy
  tests that could trigger **CPU time vs elapsed time** NOTES.

# glmbayes 0.9.4

* **Vignettes:** A vignette that previously used the `notangle` engine now
  uses the standard R Markdown vignette machinery (`knitr` /
  `rmarkdown::html_vignette`), so builds align with CRAN expectations and
  vignette index ordering should be consistent with the rest of the package.

* **OpenCL sources (`inst/cl`):** Removed unused or superseded material,
  consolidated kernels and library fragments, and aligned `.cl` layout and
  dependency tagging with the conventions used in 'openclport' and
  'nmathopencl' (prelude, shims, `nmath/` stems, family kernels under
  `src/`). See `inst/cl/README.md` for how the assembled program is stitched.

* **OpenCL program assembly:** Reworked loading so the full OpenCL program is
  built from explicit fragments (global header, `nmath` closure, family/link
  kernels) rather than ad hoc concatenation—clearer ownership of what enters
  GPU compilation and easier parity with CPU paths.

* **Tests:** Added and expanded **testthat** coverage aimed at OpenCL code
  paths (including binomial examples that exercise GPU envelope evaluation),
  complementing existing Cleveland-style checks.

* **Bug fix — binomial OpenCL:** Binomial `f2_f3` OpenCL kernels now evaluate
  the data log-likelihood with the same **proportion × trial-count**
  semantics as **`dbinom_glmb`** on the CPU (`round` successes and trials,
  clamped probability). This fixes envelope / PLSD failures for aggregated
  binomial data (e.g. `cbind(successes, failures)` / `MASS::menarche`) where
  the previous kernels treated **`y`** like a raw success count.

# glmbayes 0.9.3

* Published on CRAN.
* Version bump in response to CRAN resubmission feedback.

# glmbayes 0.9.2

* Version bump in preparation for resubmission incorporating CRAN review feedback.

# glmbayes 0.9.1

* Wrapped OpenCL-dependent examples in `\donttest{}` for CRAN compliance.
* Reduced iteration counts in rlmb Gibbs sampler example to stay within
  CRAN example time limits on slower check machines.

# glmbayes 0.9.0

First CRAN submission. This release is a stable pre-release with a
near-complete feature set relative to earlier development builds.

## Highlights

### Bayesian Generalized Linear (glmb) and Linear (lmb) modeling functions:

  `glmb()` is a Bayesian analog for the classical `glm()` function while
  `lmb()` covers Gaussian models. Calls largely mirror those for the 
  classical functions but leverage pfamilies for prior specifications.
  Method functions largely mirror those for the classical functions. 
  Samples generated by the functions are largely iid samples 
  (no MCMC convergence dignostics are needed).

### Implemented Likelihood families/ link functions:
   
  Most of the families implemented in the `glm()` function are also implemented 
  in the `glmb()` function (the `lmb()` function covers only gaussian() families). 
  Link functions that lead to log-concave likelihood functions are generally 
  implemented.  Specifically, we have the following:
  
  **Supported likelihoods:** gaussian (identity), Poisson / quasi-Poisson
  (log), binomial / quasi-binomial (logit, probit, cloglog), Gamma (log).

### Prior Family functions:

 `pfamily` constructors are used to specify priors and play the same
  kind of role for the prior specifications as `family` constructors 
  and `link` functions play for the likelihoods. Specifically, we
  have the following:

  **Supported Priors:** Normal (all families/links), Normal–Gamma and 
  independent Normal–Gamma (gaussian families), and Gamma-on-precision 
  (gaussian and Gamma families).
  
### Prior_Setup function:
 
  The package comes with a convenient `Prior_Setup()` function that provides 
  default prior input parameters for each of the implemented models. Basic calls
  (without tailoring) mirror traditional calls to the `glmb()` and `lmb()`
  functions respectively and only require the user to provide the model formula
  and (if not the gaussian family) the family/link function. 
  
  The function can also be used to easily adjust prior specifications 
  (see documentation for details).
  
### Extensive Method functions:
  
  The package comes with extensive method functions that mirror those 
  for the classical functions.  These include dedicated `print()`,
  `summary()`, `predict()` and `simulate()` functions.

### Lower Level Modeling functions:

  The package comes with lower level modeling/simulation functions
  that advanced users can use to implement block Gibbs samplers. These
  generally come with less overhead than the `glmb()` and `lmb()` functions 
  and are called internally by the the higher level modeling functions.

### RcppParallel and OpenCL GPU Acceleration Implementations
  
  Some of the simulation functions comes with use_parallel and use_opencl options
  that speed up simulation for higher dimensional models.
  
### Extensive help files, vignettes, examples and demos

  The package also comes with extensive help files for the varios functions 
  that are complemented with a rich set of vignettes. A large number of 
  examples and demos are also availabel (see the READM.md file for a sample).

---

## Earlier development history (0.1.x series)

The notes below summarize major work during the initial development series
before the 0.9.0 pre-release.

### OpenCL and GPU acceleration

- Completed the OpenCL-based grid construction framework for large models.
- Added GPU-aware envelope sizing and improved OpenCL failure handling.
- Introduced diagnostic utilities to assess OpenCL availability and
  performance.
- Improved configure scripts to detect OpenCL and provide informative
  messages.
- Expanded OpenCL documentation and added a dedicated vignette chapter.

### Parallel CPU sampling (RcppParallel)

- Enabled parallel envelope construction and parallel iid sampling.
- Added pilot functions for large-dimension grid estimation.
- Implemented thread-safe parallel sampling for independent normal-gamma
  models.

### Core statistical improvements

- Migrated to an improved independent normal-gamma simulation algorithm.
- Added theoretical derivations for independent normal-gamma regression.
- Improved UB2 and RSS minimization routines, including scaling corrections.
- Enhanced `Prior_Setup()` to support family-specific prior construction.
- Added dedicated envelope evaluation and sizing functions.

### Package infrastructure

- Significant cleanup to remove NOTES and improve CRAN readiness.
- Improved configure and Makevars files for portability.
- Added testthat tests, including OpenCL-specific tests.
- Consolidated envelope-building functions into a cleaner structure.

### Documentation

- Major updates to README and package-level documentation.
- Added multiple new vignettes and expanded existing ones.
- Improved examples for `lmb()`, `rlmb()`, and OpenCL models.

### Bug fixes (0.1.x era)

- Corrected scaling in UB2 minimization.
- Improved error handling for missing OpenCL functionality.
- Fixed various small issues uncovered during parallelization work.

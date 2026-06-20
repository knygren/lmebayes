# lmebayes (development version)

* **GLMM pilot policy in Core:** **`rglmerb()`** and **`glmerb()`** now pass
  **`gap_tol`** through to **`glmbayesCore::rGLMM()`**, which owns pilot/main
  staging (non-Gaussian default: pilot then main; Gaussian: main only).
  Removed duplicate **`n_pilot`** derivation from **`glmerb()`**.

* **Move `rLMM` to glmbayesCore:** matrix-level LMM replicate-chain
  orchestration now lives in **`glmbayesCore::rLMM()`** (v2 two-block
  driver). Re-exported from lmebayes; **`rlmerb()`** calls Core directly.

* **Move `rGLMM` to glmbayesCore:** matrix-level GLMM replicate-chain
  orchestration now lives in **`glmbayesCore::rGLMM()`** (v6 sweep-outer
  driver). Re-exported from lmebayes; **`rglmerb()`** calls Core directly.

* **Unified Block~2 return names (`fixef.*`):** `lmerb`, `glmerb`, `rlmerb`, and
  `rglmerb` now use the same `fixef.*` namespace as `rGLMM`
  (`fixef`, `fixef.mode`, `fixef.means`, `fixef.dispersion`, `fixef.iters`,
  `fixef.mu`, `fixef.init`, `pilot_chisq`).  Legacy names (`fixef_draws`,
  `coef.mode`, `tau2_draws`, `pilot_mode_test`, etc.) are removed.

* **Extract GLMM engine to `rGLMM` (glmbayesCore):** replicate-chain orchestration
  (TV calibration, pilot chi-squared, post-pilot eigenvalue upper bound,
  main-stage sweep-outer sampling) moved from `rglmerb` into
  **`glmbayesCore::rGLMM()`**, with matrix-level signature and `fixef.*`
  return layout.  `rglmerb` is now a thin `model_setup` wrapper (ICM mode,
  priors, field mapping for `glmerb`).
* **Removed legacy GLMM sampler drivers and renamed to `rglmerb`:** dropped
  development exports `rglmerb_v2`--`rglmerb_v5`, `rglmerb_experimental`, and
  internal `run_short_chains`--`run_short_chains_v5`.  The sweep-outer GLMM
  engine (formerly `rglmerb_v6`) is now `rglmerb` / class `"rglmerb"`, called
  by `glmerb()` via `glmbayesCore::run_sweep_outer_chains_v6`.

* **Two-phase Gibbs sampling in `glmerb` (new `n_pilot` argument):** for
  non-Gaussian families the ICM posterior mode is below the true posterior
  mean due to likelihood skewness (e.g. Poisson).  Starting the main sampler
  from the mode causes all stored draws to drift toward the mode, biasing the
  sample mean.  `glmerb` now runs a pilot stage of `n_pilot = 1000`
  independent short chains (default), all started at the same mode, and
  takes the column-means of their final draws as the
  estimated posterior mean (`coef.pilot.mean`), and then runs the main `n`
  draws starting from that estimate.  Because the main chain starts near the
  posterior mean from draw 1, the stored draws are near-iid in the right
  region and the sample mean is an accurate estimator of the posterior mean.
  Set `n_pilot = NULL` to restore the original single-phase behaviour.
  Ignored for `family = gaussian()` (mode = mean exactly).

* **ING `tau^2` truncation window now uses limiting-posterior quantiles
  (was: prior quantiles):** the default `disp_lower`/`disp_upper` for
  `dIndependent_Normal_Gamma` components are the 0.01/0.99 quantiles of the
  *limiting posterior* `Gamma((J+1)/2, tau2_k*(J-1)/2)` -- the weak-prior
  (`n_prior_dispersion -> 0`) limit of the Block-2 posterior Gamma for the
  precision (glmbayesCore Chapter A12, Theorem 2), inverted to the `tau^2`
  scale.  Prior-quantile windows stretch without bound as the dispersion
  prior weakens, collapsing the envelope acceptance rate while covering
  posterior mass the chain never visits; the limiting-posterior window is
  identical for every `n_prior_dispersion`, mean-matched at the classical
  `tau^2_k`, covers >= ~98% of the exact posterior for every prior
  strength (asymptotically exactly 98%), and keeps the envelope
  candidates-per-draw roughly constant as priors weaken (schools example,
  J = 47: width ratio ~2.6 vs 5.9/11.7 for the prior window at
  `pwt_dispersion` 0.2/0.1).  Derivation and trade-offs documented in
  `inst/ING_TRUNCATION_WINDOW.md`.  Note the hard truncation now clips
  ~1% of posterior mass per tail, and the higher default `disp_lower`
  makes the TV calibration less conservative (smaller `m_min`) -- still
  valid since the truncation bounds the chain's support by construction.
  Because the window no longer depends on the dispersion-prior strength,
  weak dispersion priors carry no computational penalty, so the
  `pwt_dispersion` default in `Prior_Setup_lmebayes()` remains derived
  from `pwt` (an interim flat-0.2 default, introduced earlier in this
  development cycle to keep the prior-quantile window moderate, has been
  reverted).

* **Block-2 candidate counts in `lmerb()`/`glmerb()` fits:** fits now carry
  `$iters_draws` (`n x p_re` matrix of total Block-2 candidates generated
  per stored draw, summed over the inner sweeps; from the new
  `iters_fixef_draws` output of
  `glmbayesCore::two_block_rNormal_reg_v2`) and `$iters.means` (average
  candidates per accepted draw, `colMeans(iters_draws)/m_convergence`).
  `dNormal` components are conjugate and always show 1;
  `dIndependent_Normal_Gamma` components show roughly the reciprocal
  acceptance rate of the joint `(gamma_k, tau^2_k)` envelope sampler.
  `summary()` adds a `Cand/draw` column to the Block-2 dispersion table as
  a sampler-efficiency diagnostic (values near 1-3 indicate a tight
  envelope; large values flag a mis-centered truncation window or
  ill-calibrated prior).

* **ING Gamma rate now follows the glmbayesCore default calibration
  (mean-matched):** `Prior_Setup_lmebayes()` (`ing_prior` field) and
  `pfamily_list()` previously set `rate = tau2_k * n0/2`, which matches the
  `glmbayesCore::Prior_Setup()` default `b_0 = tau2_k * (n0 + p_k - 1)/2`
  only for single-predictor components (`p_k = 1`).  For `p_k > 1` the
  implied inverse-Gamma prior on `tau^2_k` was mis-centered well below the
  classical estimate, so the default 98% truncation window
  (`disp_lower`/`disp_upper`) could exclude `hat(tau)^2_k` entirely
  (observed for a `p_k = 4` intercept component) and slow the envelope
  sampler.  The rate now uses the glmbayesCore formula; since
  `b_0 = tau2_k * (shape_ING - 1)`, the prior mean of the dispersion equals
  `tau2_k` exactly for every `n_prior_dispersion` and `p_k`, and the
  default window always brackets the classical estimate.  Covered by the
  updated `data-raw/test_pfamily_list.R` (mean-matching identity and
  window-bracketing assertions).

* **`dIndependent_Normal_Gamma` Block-2 dispersion sampling in
  `lmerb()`/`glmerb()`:** the fitters now run the pfamily-based
  `glmbayesCore::two_block_rNormal_reg_v2` sampler.  `dNormal` components
  keep the conjugate `gamma_k` draw at fixed `tau^2_k` (draws are
  bitwise-identical to the previous sampler under the same seed); ING
  components make a joint `(gamma_k, tau^2_k)` draw each inner sweep via
  the likelihood-subgradient envelope sampler, with the sampled `tau^2_k`
  fed back into the Block-1 prior precision.  The previous behavior of
  stopping after the calibration for ING components is removed; the
  conservative `disp_lower`-based TV calibration is retained.  Fits gain
  `$tau2_draws` (n x p_re matrix; constant columns for `dNormal`) and
  `$tau2.means`; `summary()` gains a per-component `tau^2` table with
  posterior mean / SD / quantiles.  Covered by
  `data-raw/test_ing_sampling.R` and the updated
  `data-raw/test_ing_calibration.R`.

* **Prior-vs-data balance guard for ING dispersion priors:** the ING
  dispersion envelope caps its log-tilt at the data contribution `J/2`
  (`J` = number of groups, the Block-2 observation count), which presumes
  a likelihood-dominated prior; a prior-dominated calibration would
  silently invalidate the envelope (biased draws were observed at small
  `J` with `pwt_dispersion` near 1).  `pfamily_list()` therefore rejects
  ING components with `n_prior_dispersion > J` (equivalently
  `pwt_dispersion > 0.5`) at construction, and the same check is enforced
  sampler-side in `glmbayesCore::two_block_rNormal_reg_v2` for hand-built
  pfamilies.  Mirrors the `n_prior <= n_w` guard added to
  `rindepNormalGamma_reg` in glmbayes/glmbayesCore.

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
  truncated support.  The fit displays the calibration
  (`conservative: ING tau^2_k = disp_lower`) and records `$convergence`
  (method `"disp_lower_bound"`, or `"<base>+disp_lower_bound"` in
  `glmerb`).  (Initially the fit stopped after the calibration; with the
  v2 sampler entry above, draws are now generated.)  On
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
  with **`rBeta_reg()`** supports Betaâ€“Binomial(identity) conjugate updates;
  **`dGamma(Inv_Dispersion = FALSE)`** with **`rGamma_Conjugate_reg()`**
  supports Gammaâ€“Poisson(identity) and Gammaâ€“Gamma(identity) rate priors.
  **`Prior_Setup()`** can calibrate conjugate hyperparameters for these
  families (weighted Poisson rate and binomial probability defaults). See
  **`?dBeta`**, **`?dGamma`**, and the Chapter 02 / Chapter 07â€“11 vignettes.

* **Vignette structure:** Reworked **Chapter 00** as a roadmap across five
  main parts plus technical appendices. **Chapter 02** is now a conceptual
  introduction to single-parameter conjugacy; worked examples move to
  **Chapter 02-S01** through **Chapter 02-S05** (Betaâ€“Binomial, Normalâ€“Normal,
  Gammaâ€“Poisson, exposure-weighted Poisson, and related topics). A **Companion
  textbooks** section in Chapter 00 indexes optional Bayes Rules! and `LearnBayes`
  appendices tied to the main GLM chapters.

* **`opencltools` import:** Core host/runtime OpenCL discovery and diagnostics
  (`detect_*`, PATH helpers, environment checks) now live in the **`opencltools`**
  package (`Imports`, >= 0.8.0). **glmbayes** keeps package-specific entry
  points (`has_opencl()`, `diagnose_glmbayes()`) that report compile-time
  OpenCL status for this build while delegating shared GPU/runtime checksâ€”reducing
  duplicated maintenance in **glmbayes**.

* **Bayes Rules! companion examples:** Optional vignette appendices reproduce
  book datasets and published posterior summaries using **`lmb()`**, **`glmb()`**,
  **`Prior_Setup()`**, and **`dNormal()`** (suggested package **`bayesrules`** for
  data only). Coverage includes **`bikes`** (Ch. 03), **`weather_perth`** (Ch. 08â€“09),
  **`equality_index`** (Ch. 10), Gammaâ€“Poisson conjugacy (Ch. 02-S04), and a
  scope note for Gamma regression (Ch. 11). Comparison tables use **printed book
  values**, not live **`rstanarm`** fits. See **Chapter 00** Â§ Companion textbooks.

* **`LearnBayes` examples:** **Chapter 02-S04**, Appendix A, maps the
  **`hearttransplants`** example from Albert (2009) / `LearnBayes` (exposure-weighted
  Gammaâ€“Poisson conjugacy) to **`glmb()`** with analytic Albert posteriors for
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
  kernels) rather than ad hoc concatenationâ€”clearer ownership of what enters
  GPU compilation and easier parity with CPU paths.

* **Tests:** Added and expanded **testthat** coverage aimed at OpenCL code
  paths (including binomial examples that exercise GPU envelope evaluation),
  complementing existing Cleveland-style checks.

* **Bug fix â€” binomial OpenCL:** Binomial `f2_f3` OpenCL kernels now evaluate
  the data log-likelihood with the same **proportion Ã— trial-count**
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

  **Supported Priors:** Normal (all families/links), Normalâ€“Gamma and 
  independent Normalâ€“Gamma (gaussian families), and Gamma-on-precision 
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

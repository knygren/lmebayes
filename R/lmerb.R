#' Bayesian linear mixed-effects model fit
#'
#' Entry point for \pkg{lmebayes} models with an \code{lmer}-like interface,
#' analogous to \code{\link{lmb}} and \code{\link{glmb}} for fixed-effects models.
#'
#' Calls \code{\link{model_setup}} on \code{formula} and \code{data} for design
#' matrices (\code{y}, \code{Z}, \code{groups}, \code{X_hyper}, etc.) and embeds
#' the resulting \code{\link[lme4]{lmer}} fit as \code{lmer}. Priors are
#' supplied as a named list of \code{\link[glmbayesCore]{pfamily}} objects
#' (\code{pfamily_list}, the Block~2 hyperpriors -- one per random-effect
#' coefficient) plus the observation-level measurement dispersion
#' (\code{dispersion_ranef}).  Both are typically built from
#' \code{\link{Prior_Setup_lmebayes}}:
#' \code{pfamily_list = pfamily_list(ps)} and either a fixed scalar
#' \code{dispersion_ranef = ps$dispersion_ranef}, a single pooled
#' \code{dGamma()} prior, or a named per-group list from
#' \code{dGamma_list(ps)}.  The Block~1 random-effect
#' covariance is reconstructed from the Block~2 pfamily dispersions
#' (\code{Sigma_ranef = diag(tau^2_k)}); \code{lmerb} does not call
#' \code{Prior_Setup_lmebayes} internally.
#'
#' Runs a two-block Gibbs sampler for \code{n} iterations. Block 1 draws
#' group-level random effects \eqn{b_j} given the current hyper means; Block 2
#' updates the hyper means (level-2 fixed effects \eqn{\boldsymbol{\gamma}_k})
#' given the current \eqn{b_j} draw, using
#' \code{\link[glmbayesCore]{multi_rNormal_reg}}
#' with the hyper design matrices from \code{design$X_hyper}.
#'
#' @details
#' \strong{Exact posterior and convergence characterisation.}
#' When variance components are treated as fixed at their \code{lmer} estimates
#' (as done here), the joint posterior over the random-effect coefficients and
#' the level-2 fixed effects is \emph{exactly} multivariate normal.  In this
#' setting the convergence of the two-block Gibbs sampler is fully
#' characterised: Corollary 1 of Nygren (2020) shows that the total variation
#' (TV) distance between the \eqn{l}-step kernel and the target density is
#' bounded above by a geometrically decreasing function of \eqn{l}, with
#' contraction rate \eqn{\lambda^* \in [0, 1)}, the maximal eigenvalue of the
#' matrix
#' \deqn{A \;=\; P_{11}^{-1/2}\,P_{12}\,P_{22}^{-1}\,P_{21}\,P_{11}^{-1/2}}
#' where \eqn{P_{11}}, \eqn{P_{22}}, and \eqn{P_{12}} are the corresponding
#' blocks of the joint precision matrix.  Because the bound is explicit and
#' computable, the required number of iterations to reach a pre-specified TV
#' tolerance \eqn{\varepsilon} can be derived analytically once \eqn{\lambda^*}
#' is known.
#'
#' \strong{TV-calibrated \code{m_convergence}.}
#' The number of inner Gibbs sweeps per stored draw (\code{m_convergence}) is
#' derived from \code{tv_tol}: \code{lmerb} computes the Remark 8 eigenvalue
#' spectrum with \code{\link[glmbayesCore]{two_block_rate_from_pfamily_list}} and inverts the
#' exact Theorem 3 bound with
#' \code{\link[glmbayesCore]{two_block_l_for_tv}}.  Because every replicate
#' chain is started at the exact joint posterior mean (computed by ICM via
#' \code{\link[glmbayesCore]{lmerb_posterior_mean}}), the mean term of the
#' bound vanishes and only the variance-convergence sum remains.  One extra
#' sweep is added because the bound applies to the block updated second in
#' each sweep (the level-2 fixed effects \eqn{\gamma}); the stored
#' random-effect draw lags by a half-step.  Each stored draw is therefore
#' guaranteed to be within \code{tv_tol} of the exact joint posterior in
#' total variation.
#'
#' @references
#' Nygren, K. (2020). \emph{On the total variation distance between multivariate
#' normal densities with applications to two-block Gibbs samplers.}
#' Unpublished manuscript.
#'
#' Jones, G. L. and Hobert, J. P. (2001). Honest exploration of intractable
#' probability distributions via Markov chain Monte Carlo.
#' \emph{Statistical Science} \bold{16}, 312--334.
#'
#' @param formula Mixed-model formula (single grouping factor; same constraints
#'   as \code{\link{model_setup}}).
#' @param data Data frame containing all variables in \code{formula}.
#' @param pfamily_list Required named list of
#'   \code{\link[glmbayesCore]{pfamily}} objects, one per random-effect
#'   coefficient (names must match the random-effect coefficient names, any
#'   order).  Supplies the Block~2 hyperpriors (\code{mu}, \code{Sigma}) and
#'   the Block~1 random-effect variances \eqn{\tau^2_k}.  \code{dNormal}
#'   components treat \eqn{\tau^2_k} (the pfamily \code{dispersion}) as
#'   known and make conjugate \eqn{\gamma_k} draws.
#'   \code{dIndependent_Normal_Gamma} components place a Gamma prior on the
#'   Block~2 precision \eqn{1/\tau^2_k}: Block~2 then makes a joint
#'   \eqn{(\gamma_k, \tau^2_k)} draw via the likelihood-subgradient envelope
#'   sampler (\code{\link[glmbayesCore]{rindepNormalGamma_reg}}), and the
#'   sampled \eqn{\tau^2_k} feeds back into the Block~1 prior precision.
#'   ING components must supply both truncation bounds: each
#'   \eqn{\tau^2_k} draw is hard-truncated to
#'   \code{[disp_lower, disp_upper]}, fixed across all inner Gibbs sweeps,
#'   with \code{disp_lower} doubling as the conservative \eqn{\tau^2_k}
#'   plug-in for the eigenvalue / TV calibration (smaller \eqn{\tau^2}
#'   increases the contraction rate \eqn{\lambda^*}, so the bound holds
#'   for every dispersion in the truncated support).  They must also
#'   satisfy the prior-vs-data guard \eqn{n_{\mathrm{prior}} \le J}
#'   (\code{pwt_dispersion} \eqn{\le 0.5}).  Typically built with
#'   \code{\link[glmbayesCore:pfamily_list.lmebayes_prior_setup]{pfamily_list}} from a
#'   \code{\link{Prior_Setup_lmebayes}} object.
#' @param dispersion_ranef Observation-level measurement dispersion
#'   \eqn{\sigma^2} for Block~1.  One of: a positive scalar (treated as
#'   known; typically \code{Prior_Setup_lmebayes(...)$dispersion_ranef}), a
#'   single \code{\link[glmbayesCore]{dGamma}()} \code{pfamily} (pooled
#'   \eqn{\sigma^2} across groups), or a named list of \code{dGamma()}
#'   objects (one per group level) from
#'   \code{\link[glmbayesCore:dGamma_list.lmebayes_prior_setup]{dGamma_list}(Prior_Setup_lmebayes(...))}.
#'   Which of these three shapes is accepted depends on \code{dispformula}
#'   (see below).
#' @param dispformula One-sided formula selecting the measurement-dispersion
#'   structure: \code{~1} (default, pooled) requires \code{dispersion_ranef}
#'   to be a fixed scalar or a single (pooled) \code{dGamma()};
#'   \code{~<group_name>}, matching the random-effects grouping factor
#'   exactly, requires \code{dispersion_ranef} to be a
#'   \code{dGamma_list(...)} (one \code{dGamma()} per group level). Any other
#'   formula is an error. \code{~1} never fits an extra reference model;
#'   \code{~<group_name>} additionally requires a \code{glmmTMB} reference fit
#'   (\pkg{glmmTMB} must be installed), stored as \code{dispersion_fit}. When
#'   \code{dispersion_ranef = dGamma_list(Prior_Setup_lmebayes(..., dispformula
#'   = dispformula))}, that call already fit this reference model to
#'   calibrate the priors, and \code{lmerb()} reuses it here rather than
#'   fitting \code{glmmTMB} a second time; keep \code{dispformula} identical
#'   between the two calls, since it is not re-validated against the reused
#'   fit. \code{lmer} is always the plain \code{\link[lme4]{lmer}} fit
#'   regardless of \code{dispformula}; the sampler route (pooled vs.
#'   per-group) already follows from \code{dispersion_ranef}'s shape alone.
#' @param n Number of iid draws per group (default \code{1000L}, as in \code{\link{lmb}}).
#' @param tv_tol Total variation tolerance per stored draw, in (0, 1)
#'   (default \code{0.01}, the conventional threshold of the honest-burn-in
#'   literature; Jones and Hobert 2001).  The number of inner Gibbs sweeps
#'   per stored draw is derived so that each draw is within \code{tv_tol} of
#'   the exact joint posterior in total variation (Nygren 2020, Theorem 3;
#'   see Details).  To certify the whole \code{n}-draw sample at level
#'   \eqn{\alpha} pass \code{tv_tol = alpha / n}; the cost grows only
#'   logarithmically in \code{1/tv_tol}.
#' @param gap_tol Legacy mode--mean gap tolerance for the pilot stage when
#'   any Block~2 component uses \code{dIndependent_Normal_Gamma} and
#'   \code{tv_tol} is \code{NULL}.  When \code{tv_tol} is set (default),
#'   the pilot chain count is chosen by cost optimization instead.  Set
#'   \code{NULL} to skip the pilot unless \code{tv_tol} is set.  Ignored for
#'   all-\code{dNormal} models.
#' @param mode_gap_max Maximum per-coordinate mode--mean gap (in posterior
#'   standard deviation units) used to calibrate pilot inner sweeps when ING
#'   Block~2 components run a pilot stage (default \code{1.0}).  Ignored for
#'   all-\code{dNormal} models.
#' @param REML Logical; passed to \code{\link{model_setup}}.
#' @param control \code{\link[lme4]{lmerControl}} settings; passed to \code{model_setup}.
#' @param start Optional starting values; passed to \code{model_setup}.
#' @param verbose Verbosity flag; passed to \code{model_setup}.
#' @param subset Optional subset; passed to \code{model_setup}.
#' @param weights Optional weights; passed to \code{model_setup}.
#' @param na.action Missing-data handler; passed to \code{model_setup}.
#' @param offset Optional offset; passed to \code{model_setup}.
#' @param contrasts Optional contrasts; passed to \code{model_setup}.
#' @param devFunOnly If \code{TRUE}, return deviance function only; passed to \code{model_setup}.
#' @param simulate Logical (default \code{TRUE}).  When \code{TRUE} the
#'   two-block Gibbs sampler is run for \code{n} iterations and posterior draws
#'   are stored.  When \code{FALSE} only the ICM algorithm is run: the exact
#'   posterior means (\code{fixef.mode}, \code{ranef.mode}) are computed and
#'   returned immediately without any sampling.  Simulation-only fields
#'   (\code{coefficients}, \code{fixef.means}, \code{fixef}) are
#'   \code{NULL} when \code{simulate = FALSE}.
#' @param fixef Optional named list of hyper-parameter vectors (Block 2 state).
#'   When \code{NULL} (default), iter-0 means are taken from the
#'   \code{pfamily_list} prior means.
#' @param diag_sweeps Temporary diagnostic flag (ING models with pilot).
#'   Non-\code{dNormal} sampling already runs via
#'   \code{glmbayesCore::run_sweep_outer_chains_v6()} (R sweep-outer;
#'   pilot then main).  When \code{TRUE}, each stage auto-prints one combined
#'   Block~2 chain-mean table when that stage finishes (same layout as
#'   \code{print()} on \code{$sweep_history}).  \code{sweep_history} is
#'   collected regardless.  Default \code{FALSE}.
#' @param progbar Logical. Show text progress bars during sampling (passed to
#'   \code{\link{rlmerb}}). Default \code{NULL}: \code{FALSE} when
#'   \code{diag_sweeps = TRUE}, otherwise \code{TRUE}.
#' @param ... Reserved for future use.
#' @return Object of class \code{"lmerb"}: a list with the following
#'   components (parallel to \code{\link{glmb}} and \code{\link{lmb}}):
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{formula}}{The formula supplied.}
#'     \item{\code{lmer}}{\code{\link[lme4]{lmer}} fit from
#'       \code{model_setup} (full \code{formula}), embedded as a sub-object —
#'       analogous to \code{glmb$glm} and \code{lmb$lm}.  Use
#'       \code{coef(fit$lmer)} for per-group classical coefficients.  Always
#'       the plain pooled-dispersion \code{lmer} fit, regardless of
#'       \code{dispformula}.}
#'     \item{\code{dispformula}}{The \code{dispformula} supplied.}
#'     \item{\code{dispersion_fit}}{\code{NULL} when \code{dispformula = ~1};
#'       otherwise the diagnostic-only \code{\link[glmmTMB]{glmmTMB}} fit with
#'       per-group residual dispersion (\code{dispformula}).}
#'     \item{\code{prior}}{Normalized prior container: \code{pfamily_list}
#'       (as supplied, reordered to the RE coefficient names),
#'       \code{dispersion_ranef}, the reconstructed \code{Sigma_ranef}, and
#'       the per-component \code{prior_list} (\code{mu_fixef},
#'       \code{Sigma_fixef}, \code{dispersion_fixef}) — analogous to
#'       \code{glmb$Prior}.}
#'     \item{\code{model_setup}}{The \code{\link{model_setup}} object built
#'       inside \code{lmerb} from \code{formula} and \code{data}.}
#'     \item{\code{fixef.mode}}{Named list of exact posterior mode (= mean,
#'       since the joint posterior is Gaussian) vectors for the level-2 fixed
#'       effects \eqn{\gamma_k}, computed by
#'       \code{\link[glmbayesCore]{lmerb_posterior_mean}} (ICM).}
#'     \item{\code{ranef.mode}}{\eqn{J \times p_{\mathrm{re}}} numeric matrix
#'       of exact posterior mode random effects from ICM.  Rows are group
#'       levels (\code{levels(design$groups)}); columns are
#'       \code{design$re_coef_names}.}
#'     \item{\code{fixef.means}}{Named list of posterior mean vectors computed
#'       as \code{colMeans(fixef[[k]])} — the MCMC estimate of the
#'       level-2 fixed effects.  \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{fixef}}{Named list of \eqn{n \times q_k} matrices of
#'       Block 2 draws, one per RE component.  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{coefficients}}{\code{data.frame} with \code{n * J} rows:
#'       \code{draw}, the grouping-factor column, and one column per RE
#'       variable.  Average over \code{draw} within each group for posterior
#'       means (see Examples).  \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{fixef.dispersion}}{\eqn{n \times p_{\mathrm{re}}} matrix of the
#'       Block~2 dispersion (\eqn{\tau^2_k}) at each stored draw: sampled
#'       values for \code{dIndependent_Normal_Gamma} components, constant
#'       columns (the fixed \code{dispersion}) for \code{dNormal} components.
#'       \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{fixef.dispersion.mean}}{Named vector of posterior means of
#'       \eqn{\tau^2_k} (\code{colMeans(fixef.dispersion)}).  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{sigma2}}{Observation-level residual variance \eqn{\sigma^2}:
#'       a scalar when \code{dispersion_ranef} is fixed, a length-\code{n} vector
#'       of final-sweep draws when \code{dispersion_ranef} is \code{dGamma()};
#'       \code{NULL} when not applicable.  When \code{simulate = FALSE}, the
#'       prior plug-in scalar (\code{prior$dispersion_ranef}) is returned for
#'       Gaussian models.}
#'     \item{\code{sigma2.mean}}{Posterior mean of \eqn{\sigma^2}; equals
#'       \code{sigma2} when fixed.  \code{NULL} when \code{sigma2} is
#'       \code{NULL}.}
#'     \item{\code{draw_engine}}{Name of the Block~1 sampling engine used
#'       (e.g. \code{"rGLMM_sweep_ing_block1_ind"}). \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{fixef.iters}}{\eqn{n \times p_{\mathrm{re}}} matrix of the
#'       total number of Block~2 candidates generated per stored draw,
#'       summed over the \code{m_convergence} inner sweeps.
#'       \code{dIndependent_Normal_Gamma} components count envelope
#'       accept-reject candidates until acceptance (as \code{iters} in
#'       \pkg{glmbayes} samplers); \code{dNormal} components count exactly
#'       one conjugate draw per sweep.  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{fixef.iters.mean}}{Named vector of average candidates per
#'       accepted Block~2 draw
#'       (\code{colMeans(fixef.iters)/m_convergence}; equals 1 for
#'       \code{dNormal} components, and approximately the reciprocal
#'       acceptance rate of the ING envelope sampler otherwise).
#'       \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{fixef.mu}}{Numeric matrix \code{p_re x J} of Block 1 prior
#'       means at the final Gibbs state (from
#'       \code{\link[glmbayesCore]{build_mu_all}}).}
#'     \item{\code{convergence}}{List describing the sweep-count calibration:
#'       \code{method} (\code{"exact"}, or \code{"local_gaussian_mode"} for
#'       non-Gaussian \code{\link{glmerb}}), \code{tv_tol},
#'       \code{lambda_star}, \code{eigenvalues}, \code{m_min} (derived
#'       minimum sweeps), and \code{m_convergence} (sweeps actually used).
#'       \code{NULL} when \code{simulate = FALSE}.}
#'   }
#' @example inst/examples/Ex_lmerb.R
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{model_setup}},
#'   \code{\link[glmbayesCore]{build_mu_all}},
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg}},
#'   \code{\link[glmbayesCore]{lmerb_posterior_mean}},
#'   \code{\link[glmbayesCore]{block_rNormalReg}},
#'   \code{\link{lmb}}, \code{\link{glmb}}
#' @param digits Number of significant digits to use when printing.
#' @title Fit a Bayesian linear mixed-effects model (LMM) to data, via two-Block Gibbs sampling
#' @aliases lmerb print.lmerb
#' @export
lmerb <- function(
    formula,
    data = NULL,
    pfamily_list,
    dispersion_ranef,
    dispformula = ~1,
    n = 1000L,
    tv_tol = 0.01,
    gap_tol = 0.0196,
    mode_gap_max = 1.0,
    diag_sweeps = FALSE,
    progbar = NULL,
    simulate = TRUE,
    REML = TRUE,
    control = lme4::lmerControl(),
    start = NULL,
    verbose = 0L,
    subset,
    weights,
    na.action,
    offset,
    contrasts = NULL,
    devFunOnly = FALSE,
    fixef = NULL,
    ...
) {
  cl <- match.call()
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (missing(pfamily_list) || is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required. Build it with ",
      "pfamily_list(Prior_Setup_lmebayes(...)) and pass the result to lmerb().",
      call. = FALSE
    )
  }
  if (missing(dispersion_ranef)) {
    stop(
      "'dispersion_ranef' is required for lmerb(). Typically ",
      "Prior_Setup_lmebayes(...)$dispersion_ranef.",
      call. = FALSE
    )
  }

  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) {
    stop("'n' must be at least 1.", call. = FALSE)
  }
  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }
  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
  }

  setup_args <- list(
    formula = formula,
    data = data,
    REML = REML,
    control = control,
    verbose = verbose,
    devFunOnly = devFunOnly
  )
  if (!missing(start) && !is.null(start)) {
    setup_args$start <- start
  }
  if (!missing(subset)) {
    setup_args$subset <- subset
  }
  if (!missing(weights)) {
    setup_args$weights <- weights
  }
  if (!missing(na.action)) {
    setup_args$na.action <- na.action
  }
  if (!missing(offset)) {
    setup_args$offset <- offset
  }
  if (!missing(contrasts)) {
    setup_args$contrasts <- contrasts
  }

  design <- do.call(model_setup, c(setup_args, list(...)))
  if (!inherits(design, "model_setup")) {
    stop("model_setup() must return a model_setup object.", call. = FALSE)
  }

  prior <- glmbayesCore:::.lmebayes_priors_from_pfamily_list(
    pfamily_list     = pfamily_list,
    dispersion_ranef = dispersion_ranef,
    design           = design,
    family           = gaussian(),
    fn_name          = "lmerb"
  )

  dispformula_kind <- .lmebayes_validate_dispformula(
    dispformula = dispformula,
    group_name  = design$group_name,
    family      = gaussian(),
    disp_mode   = prior$dispersion_mode
  )
  dispersion_fit <- NULL
  if (identical(dispformula_kind, "group")) {
    ## dGamma_list(Prior_Setup_lmebayes(..., dispformula = dispformula))
    ## already carries its glmmTMB reference fit forward as an attribute;
    ## reuse it instead of re-fitting glmmTMB here.
    dispersion_fit <- attr(dispersion_ranef, "dispersion_fit")
    if (is.null(dispersion_fit)) {
      dispersion_fit <- .lmebayes_fit_glmmtmb_dispersion(
        formula           = formula,
        data              = data,
        family            = gaussian(),
        dispformula       = dispformula,
        REML              = REML,
        mer_optional_args = glmbayesCore:::.lmebayes_mer_optional_args(
          start     = start,
          subset    = subset,
          weights   = weights,
          na.action = na.action,
          offset    = offset,
          contrasts = contrasts
        )
      )
    }
  }

  lmer_fit <- design$lmer_fit

  if (is.null(fixef)) {
    fixef <- lapply(prior$prior_list, `[[`, "mu_fixef")
    names(fixef) <- design$re_coef_names
  }

  if (!isTRUE(simulate)) {
    icm <- .lmebayes_icm_at_fixed_vc(
      design = design,
      prior  = prior,
      family = gaussian()
    )
    .lmebayes_print_icm_simulate_false(
      prior    = prior,
      re_names = design$re_coef_names,
      icm      = icm,
      header   = "--- lmerb: Block 2 fixed effects ---"
    )
    return(structure(
      list(
        call         = cl,
        formula      = formula,
        lmer         = lmer_fit,
        dispformula    = dispformula,
        dispersion_fit = dispersion_fit,
        prior        = prior,
        model_setup  = design,
        fixef.mode   = icm$fixef,
        fixef.init   = icm$fixef_init,
        ranef.mode   = icm$b_mean,
        fixef.means  = NULL,
        fixef        = NULL,
        coefficients = NULL,
        joint_mode   = icm$joint_mode,
        sigma2       = icm$sigma2,
        sigma2.mean  = icm$sigma2,
        sigma2.mode  = icm$sigma2,
        tau2.mode    = icm$tau2,
        fixef.mu     = as.matrix(
          glmbayesCore::build_mu_all(design, icm$fixef)$mu_all
        ),
        draw_engine  = NULL
      ),
      class = c("lmerb", "list")
    ))
  }

  # ICM posterior mean, block1_prior, convergence calibration, and sampling
  # are all handled inside rlmerb.
  sampler <- glmbayesCore::rlmerb(
    n               = n,
    design          = design,
    prior           = prior,
    dispersion_ranef = dispersion_ranef,
    tv_tol        = tv_tol,
    progbar       = progbar,
    verbose       = TRUE,
    gap_tol             = gap_tol,
    mode_gap_max        = mode_gap_max,
    diag_sweeps         = diag_sweeps
  )

  convergence_info <- sampler$convergence
  m_convergence    <- sampler$m_convergence
  run_pilot <- !is.null(sampler$n_pilot) && sampler$n_pilot > 0L

  structure(
    list(
      call                  = cl,
      formula               = formula,
      lmer                  = lmer_fit,
      dispformula           = dispformula,
      dispersion_fit        = dispersion_fit,
      prior                 = prior,
      model_setup           = design,
      fixef.mode            = sampler$fixef.mode,
      fixef.init            = if (run_pilot) sampler$fixef.init else NULL,
      ranef.mode            = sampler$ranef.mode,
      fixef.means           = sampler$fixef.means,
      fixef                 = sampler$fixef,
      coefficients          = sampler$coefficients,
      fixef.dispersion      = sampler$fixef.dispersion,
      fixef.dispersion.mean = sampler$fixef.dispersion.mean,
      fixef.iters           = sampler$fixef.iters,
      fixef.iters.mean      = sampler$fixef.iters.mean,
      ranef.iters           = sampler$ranef.iters,
      ranef.iters.mean      = sampler$ranef.iters.mean,
      sigma2                = sampler$sigma2,
      sigma2.mean           = sampler$sigma2.mean,
      sigma2.iters          = sampler$sigma2.iters,
      sigma2.iters.mean     = sampler$sigma2.iters.mean,
      fixef.mu              = sampler$fixef.mu,
      draw_engine           = sampler$draw_engine,
      m_convergence         = m_convergence,
      pilot_chisq           = sampler$pilot_chisq,
      gap_tol               = if (isTRUE(prior$any_non_normal)) gap_tol else NULL,
      mode_gap_max          = if (isTRUE(prior$any_non_normal)) mode_gap_max else NULL,
      convergence           = convergence_info,
      sweep_history         = list(
        pilot = if (run_pilot && !is.null(sampler$pilot$sweep_history)) {
          sampler$pilot$sweep_history
        } else {
          NULL
        },
        main = sampler$sweep_history
      ),
      pilot                 = sampler$pilot
    ),
    class = c("lmerb", "list")
  )
}

#' @rdname lmerb
#' @method print lmerb
#' @param x Object of class \code{"lmerb"}.
#' @param sweep_history If \code{TRUE}, print stored Block~2 sweep history
#'   (see \code{$sweep_history}) after the usual summary: one combined table
#'   per stage with mode plus all inner sweeps.
#' @param sweep_history_stage Which stage to print when \code{sweep_history =
#'   TRUE}: \code{"main"}, \code{"pilot"}, or \code{"both"}.
#' @param max_sweeps Passed to \code{print()} on sweep-history objects; limits
#'   how many inner sweeps are shown per stage.
#' @param components Optional RE components passed to sweep-history \code{print()}.
#' @param covariate Optional covariate names passed to sweep-history \code{print()}.
#' @export
print.lmerb <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    sweep_history = FALSE,
    sweep_history_stage = c("main", "pilot", "both"),
    max_sweeps = Inf,
    components = NULL,
    covariate = NULL,
    ...
) {

  re_names <- x$model_setup$re_coef_names
  grp      <- x$model_setup$group_name
  n_obs    <- length(x$model_setup$y)
  n_grp    <- nlevels(x$model_setup$groups)
  simulated <- !is.null(x$coefficients)

  # --- Call ---
  cat("Call:\n  ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  # --- Header line ---
  if (simulated) {
    n_draws <- nrow(x$fixef[[re_names[1L]]])
    cat(sprintf(
      "Bayesian linear mixed model  [%d draws, two-block Gibbs]\n", n_draws))
  } else {
    cat("Bayesian linear mixed model  [ICM only; use simulate = TRUE for draws]\n")
  }
  cat("Formula:", deparse1(x$formula), "\n\n")

  # --- Variance components ---
  any_non_normal <- isTRUE(x$any_non_normal) || isTRUE(x$prior$any_non_normal)
  if (any_non_normal) {
    cat("Random effects (lmer reference; tau^2 sampled for non-dNormal components):\n")
  } else {
    cat("Random effects (variance components fixed at lmer estimates):\n")
  }
  print(lme4::VarCorr(x$lmer), comp = "Std.Dev.", digits = digits)
  cat(sprintf("Number of obs: %d,  groups: %s, %d\n\n", n_obs, grp, n_grp))
  if (any_non_normal && !is.null(x$fixef.dispersion.mean)) {
    cat("Posterior mean tau^2_k: ",
        paste(sprintf("%s = %.4g", names(x$fixef.dispersion.mean),
                      x$fixef.dispersion.mean),
              collapse = ", "),
        "\n\n", sep = "")
  }

  # --- Block 2 hyperparameter reference (ICM) vs MCMC means when simulated ---
  mode_col <- if (any_non_normal) {
    cat("--- Block 2 hyperparameters (gamma at lmer tau^2 plug-in; MCMC means when simulated) ---\n\n")
    "gamma @ lmer tau2"
  } else {
    cat("--- Posterior means (ICM exact, under fixed variance components) ---\n\n")
    "fixef.mode"
  }

  rows <- do.call(rbind, lapply(re_names, function(k) {
    nms <- names(x$fixef.mode[[k]])
    data.frame(
      re  = k,
      par = nms,
      mode = unname(x$fixef.mode[[k]]),
      stringsAsFactors = FALSE
    )
  }))

  w_re  <- max(nchar(rows$re),  nchar("RE component"))
  w_par <- max(nchar(rows$par), nchar("parameter"))

  if (!simulated) {
    cat(sprintf("  %-*s  %-*s  %12s\n",
                w_re, "RE component", w_par, "parameter", mode_col))
    cat(sprintf("  %s  %s  %s\n",
                strrep("-", w_re), strrep("-", w_par), strrep("-", 12L)))
    for (i in seq_len(nrow(rows))) {
      cat(sprintf("  %-*s  %-*s  %12.*f\n",
                  w_re, rows$re[i], w_par, rows$par[i],
                  digits, rows$mode[i]))
    }
    cat("\n")

  } else {
    rows$means <- unlist(lapply(re_names, function(k) unname(x$fixef.means[[k]])))
    rows$sd    <- unlist(lapply(re_names, function(k) {
      apply(x$fixef[[k]], 2L, sd)
    }))

    cat(sprintf("  %-*s  %-*s  %12s  %12s  %10s\n",
                w_re, "RE component", w_par, "parameter",
                mode_col, "fixef.means", "draws SD"))
    cat(sprintf("  %s  %s  %s  %s  %s\n",
                strrep("-", w_re), strrep("-", w_par),
                strrep("-", 12L), strrep("-", 12L), strrep("-", 10L)))
    for (i in seq_len(nrow(rows))) {
      cat(sprintf("  %-*s  %-*s  %12.*f  %12.*f  %10.*f\n",
                  w_re, rows$re[i], w_par, rows$par[i],
                  digits, rows$mode[i],
                  digits, rows$means[i],
                  digits, rows$sd[i]))
    }
    cat("\n")
  }

  if (isTRUE(sweep_history) && !is.null(x$sweep_history)) {
    stage <- match.arg(sweep_history_stage)
    stages <- switch(
      stage,
      main = "main",
      pilot = "pilot",
      both = c("pilot", "main")
    )
    for (st in stages) {
      hist <- x$sweep_history[[st]]
      if (!is.null(hist)) {
        print(
          hist,
          max_sweeps = max_sweeps,
          components = components,
          covariate  = covariate,
          digits     = digits,
          ...
        )
      }
    }
  }

  invisible(x)
}

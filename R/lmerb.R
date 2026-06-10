#' Bayesian linear mixed model fit (draft)
#'
#' Entry point for \pkg{lmebayes} models with an \code{lmer}-like interface,
#' analogous to \code{\link{lmb}} and \code{\link{glmb}} for fixed-effects models.
#'
#' Calls \code{\link{model_setup}} on \code{formula} and \code{data} for design
#' matrices (\code{y}, \code{Z}, \code{groups}, \code{X_hyper}, etc.) and embeds
#' the resulting \code{\link[lme4]{lmer}} fit as \code{lmer}. Measurement priors
#' (\code{dispersion_ranef}, \code{Sigma_ranef}, Block~2 \code{prior_list}) must
#' be supplied via \code{measurement_prior_list} from a prior call to
#' \code{\link{Prior_Setup_lmebayes}}; \code{lmerb} does not call
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
#' spectrum with \code{\link[glmbayesCore]{two_block_rate}} and inverts the
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
#'   as \code{\link{model_setup}}). Must match
#'   \code{measurement_prior_list$formula}.
#' @param data Data frame containing all variables in \code{formula}.
#' @param measurement_prior_list Required object from
#'   \code{\link{Prior_Setup_lmebayes}}. Supplies \code{dispersion_ranef},
#'   \code{Sigma_ranef}, and Block~2 priors via \code{prior_list} (including
#'   iter-0 \code{mu_fixef} means). Call \code{\link{Prior_Setup_lmebayes}}
#'   explicitly before \code{lmerb}; design matrices come from
#'   \code{\link{model_setup}} inside \code{lmerb}, not from this object.
#' @param n Number of iid draws per group (default \code{1000L}, as in \code{\link{lmb}}).
#' @param tv_tol Total variation tolerance per stored draw, in (0, 1)
#'   (default \code{0.01}, the conventional threshold of the honest-burn-in
#'   literature; Jones and Hobert 2001).  The number of inner Gibbs sweeps
#'   per stored draw is derived so that each draw is within \code{tv_tol} of
#'   the exact joint posterior in total variation (Nygren 2020, Theorem 3;
#'   see Details).  To certify the whole \code{n}-draw sample at level
#'   \eqn{\alpha} pass \code{tv_tol = alpha / n}; the cost grows only
#'   logarithmically in \code{1/tv_tol}.
#' @param m_convergence Optional integer override for the number of inner
#'   Gibbs sweeps per stored draw.  When \code{NULL} (default) the
#'   \code{tv_tol}-derived value is used.  A supplied value is floored at the
#'   derived minimum: \code{max(m_convergence, m_min)} is used, with a
#'   warning if the value had to be raised.
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
#'   posterior means (\code{coef.mode}, \code{ranef.mode}) are computed and
#'   returned immediately without any sampling.  Simulation-only fields
#'   (\code{coefficients}, \code{coef.means}, \code{fixef_draws}) are
#'   \code{NULL} when \code{simulate = FALSE}.
#' @param fixef Optional named list of hyper-parameter vectors (Block 2 state).
#'   When \code{NULL} (default), iter-0 means are taken from
#'   \code{measurement_prior_list$prior_list$mu_fixef}.
#' @param seed Optional; sets the RNG seed before sampling.
#' @param ... Reserved for future use.
#' @return Object of class \code{"lmerb"}: a list with the following
#'   components (parallel to \code{\link{glmb}} and \code{\link{lmb}}):
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{formula}}{The formula supplied.}
#'     \item{\code{lmer}}{\code{\link[lme4]{lmer}} fit from
#'       \code{model_setup} (full \code{formula}), embedded as a sub-object —
#'       analogous to \code{glmb$glm} and \code{lmb$lm}.  Use
#'       \code{coef(fit$lmer)} for per-group classical coefficients.}
#'     \item{\code{prior}}{The \code{measurement_prior_list} object supplied by
#'       the caller (from \code{\link{Prior_Setup_lmebayes}}), stored for
#'       reference and re-use — analogous to \code{glmb$Prior}.}
#'     \item{\code{model_setup}}{The \code{\link{model_setup}} object built
#'       inside \code{lmerb} from \code{formula} and \code{data}.}
#'     \item{\code{coef.mode}}{Named list of exact posterior mode (= mean,
#'       since the joint posterior is Gaussian) vectors for the level-2 fixed
#'       effects \eqn{\gamma_k}, computed by
#'       \code{\link[glmbayesCore]{lmerb_posterior_mean}}
#'       (ICM).  Analogous to \code{glmb$coef.mode}.}
#'     \item{\code{ranef.mode}}{\eqn{J \times p_{\mathrm{re}}} numeric matrix
#'       of exact posterior mode random effects from ICM.  Rows are group
#'       levels (\code{levels(design$groups)}); columns are
#'       \code{design$re_coef_names}.}
#'     \item{\code{coef.means}}{Named list of posterior mean vectors computed
#'       as \code{colMeans(fixef_draws[[k]])} — the MCMC estimate of the
#'       level-2 fixed effects.  Analogous to \code{glmb$coef.means}.
#'       \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{fixef_draws}}{Named list of \eqn{n \times q_k} matrices of
#'       Block 2 draws, one per RE component.  \code{NULL} when
#'       \code{simulate = FALSE}.}
#'     \item{\code{coefficients}}{\code{data.frame} with \code{n * J} rows:
#'       \code{draw}, the grouping-factor column, and one column per RE
#'       variable.  Average over \code{draw} within each group for posterior
#'       means (see Examples).  \code{NULL} when \code{simulate = FALSE}.}
#'     \item{\code{mu_all}}{Numeric matrix \code{p_re x J} of Block 1 prior
#'       means at the final Gibbs state (from
#'       \code{\link[glmbayesCore]{build_mu_all}}).}
#'     \item{\code{convergence}}{List describing the sweep-count calibration:
#'       \code{method} (\code{"exact"}, or \code{"local_gaussian_mode"} for
#'       non-Gaussian \code{\link{glmerb}}), \code{tv_tol},
#'       \code{lambda_star}, \code{eigenvalues}, \code{m_min} (derived
#'       minimum sweeps), and \code{m_convergence} (sweeps actually used).
#'       \code{NULL} when \code{simulate = FALSE}.}
#'   }
#' @examples
#' \donttest{
#'   source(system.file("examples", "Ex_lmerb.R", package = "lmebayes"))
#' }
#' @seealso \code{\link{Prior_Setup_lmebayes}}, \code{\link{model_setup}},
#'   \code{\link[glmbayesCore]{build_mu_all}},
#'   \code{\link[glmbayesCore]{two_block_rNormal_reg}},
#'   \code{\link[glmbayesCore]{lmerb_posterior_mean}},
#'   \code{\link[glmbayesCore]{block_rNormalReg}},
#'   \code{\link{lmb}}, \code{\link{glmb}}
#' @export
lmerb <- function(
    formula,
    data = NULL,
    measurement_prior_list,
    n = 1000L,
    tv_tol = 0.01,
    m_convergence = NULL,
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
    seed = NULL,
    ...
) {
  cl <- match.call()
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (is.null(data) || !is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (missing(measurement_prior_list) || is.null(measurement_prior_list)) {
    stop(
      "'measurement_prior_list' is required. ",
      "Build it with Prior_Setup_lmebayes() and pass the result to lmerb().",
      call. = FALSE
    )
  }
  if (!inherits(measurement_prior_list, "lmebayes_prior_setup")) {
    stop(
      "'measurement_prior_list' must be an object from Prior_Setup_lmebayes().",
      call. = FALSE
    )
  }
  if (!identical(deparse(formula), deparse(measurement_prior_list$formula))) {
    stop(
      "'formula' does not match measurement_prior_list$formula.",
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
  if (!is.null(m_convergence)) {
    if (!is.numeric(m_convergence) || length(m_convergence) != 1L ||
        !is.finite(m_convergence) || m_convergence < 1) {
      stop("'m_convergence' must be NULL or a single integer >= 1.",
           call. = FALSE)
    }
    m_convergence <- as.integer(m_convergence)
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

  if (!identical(design$re_coef_names, names(measurement_prior_list$prior_list))) {
    stop(
      "measurement_prior_list$prior_list names must match design$re_coef_names.",
      call. = FALSE
    )
  }

  lmer_fit <- design$lmer_fit

  if (is.null(fixef)) {
    fixef <- lapply(measurement_prior_list$prior_list, `[[`, "mu_fixef")
    names(fixef) <- design$re_coef_names
  }

  # Common starting state for every inner Gibbs run.  The TV bound in Nygren
  # (2020) Corollary 1 is tightest when the chain starts at the joint posterior
  # mean; using that point minimises the epsilon achievable in m_convergence
  # steps.  Use lmerb_posterior_mean() to find the exact posterior mean via ICM.
  fixef_lmer <- fixef   # lmer-derived starting values (for diagnostic printing)
  pm <- glmbayesCore::lmerb_posterior_mean(design, measurement_prior_list)
  fixef_start <- pm$fixef

  # Diagnostic table: lmer start vs ICM posterior mean, one row per parameter
  hdr <- sprintf("  %-18s  %-30s  %12s  %12s",
                 "RE component", "parameter", "lmer (start)", "post mean (ICM)")
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat("--- lmerb: Block 2 fixed effects ---\n")
  cat(hdr, "\n")
  cat(sep, "\n")
  for (k in design$re_coef_names) {
    nms_k  <- names(fixef_lmer[[k]])
    lmer_v <- fixef_lmer[[k]]
    pm_v   <- fixef_start[[k]]
    for (nm in nms_k) {
      cat(sprintf("  %-18s  %-30s  %12.4f  %12.4f\n",
                  k, nm, lmer_v[[nm]], pm_v[[nm]]))
    }
  }
  cat(sprintf("  (ICM converged: %s, %d iter, delta = %.2e)\n\n",
              pm$converged, pm$iterations, pm$delta))

  if (is.null(measurement_prior_list$dispersion_ranef)) {
    stop(
      "'measurement_prior_list' must contain 'dispersion_ranef' for lmerb().",
      call. = FALSE
    )
  }
  block1_prior <- .lmebayes_block1_prior_list(measurement_prior_list)

  # When simulate=FALSE return only the ICM posterior means immediately.
  if (!isTRUE(simulate)) {
    return(structure(
      list(
        call        = cl,
        formula     = formula,
        lmer        = lmer_fit,
        prior       = measurement_prior_list,
        model_setup = design,
        coef.mode   = fixef_start,
        ranef.mode  = pm$b_mean,
        coef.means  = NULL,
        fixef_draws = NULL,
        coefficients = NULL,
        mu_all      = as.matrix(
          glmbayesCore::build_mu_all(design, fixef_start)$mu_all
        )
      ),
      class = c("lmerb", "list")
    ))
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

  # Block 2 prior: fixed hyperpriors from Prior_Setup_lmebayes.
  block2_prior_list <- stats::setNames(
    lapply(re_names, function(k) {
      pl_k <- measurement_prior_list$prior_list[[k]]
      list(
        mu         = pl_k$mu_fixef,
        Sigma      = pl_k$Sigma_fixef,
        dispersion = pl_k$dispersion_fixef
      )
    }),
    re_names
  )

  # TV-calibrated number of inner Gibbs sweeps per stored draw.  With fixed
  # variance components the joint posterior is exactly multivariate normal,
  # so the Remark 8 spectrum (Nygren 2020) gives the exact Theorem 3 TV bound
  # for the l-step kernel.  Every replicate chain starts at the joint
  # posterior mean (fixef_start, via ICM), so the mean term vanishes (D0 = 0).
  # The bound applies to the block updated second (gamma); the stored b draw
  # lags by a half-step, hence the + 1L.
  rate <- glmbayesCore::two_block_rate(
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = block1_prior,
    prior_list_block2 = block2_prior_list,
    family            = gaussian(),
    group_levels      = group_levels
  )
  m_min <- glmbayesCore::two_block_l_for_tv(
    rate, tv_tol, method = "theorem3"
  ) + 1L
  if (is.null(m_convergence)) {
    m_convergence <- m_min
  } else if (m_convergence < m_min) {
    warning(
      "lmerb: m_convergence = ", m_convergence, " is below the derived ",
      "minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
      "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }
  cat(sprintf(
    "--- lmerb: convergence calibration [exact]: lambda* = %.4f, tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n",
    rate$lambda_star, tv_tol, m_min, m_convergence
  ))
  convergence_info <- list(
    method        = "exact",
    tv_tol        = tv_tol,
    lambda_star   = rate$lambda_star,
    eigenvalues   = rate$eigenvalues,
    m_min         = m_min,
    m_convergence = m_convergence
  )

  sampler <- glmbayesCore::two_block_rNormal_reg(
    n                 = n,
    y                 = design$y,
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = block1_prior,
    prior_list_block2 = block2_prior_list,
    fixef_start       = fixef_start,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = design$group_name,
    family            = gaussian(),
    m_convergence     = m_convergence,
    seed              = seed,
    progbar           = TRUE
  )

  structure(
    list(
      call         = cl,
      formula      = formula,
      lmer         = lmer_fit,
      prior        = measurement_prior_list,
      model_setup  = design,
      coef.mode    = fixef_start,
      ranef.mode   = pm$b_mean,
      coef.means   = lapply(sampler$fixef_draws, colMeans),
      fixef_draws  = sampler$fixef_draws,
      coefficients = sampler$coefficients,
      mu_all       = sampler$mu_all_last,
      convergence  = convergence_info
    ),
    class = c("lmerb", "list")
  )
}

#' Print method for lmerb objects
#'
#' @param x Object of class \code{"lmerb"}.
#' @param digits Number of significant digits (default
#'   \code{max(3, getOption("digits") - 3)}).
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.lmerb <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

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
    n_draws <- nrow(x$fixef_draws[[re_names[1L]]])
    cat(sprintf(
      "Bayesian linear mixed model  [%d draws, two-block Gibbs]\n", n_draws))
  } else {
    cat("Bayesian linear mixed model  [ICM only; use simulate = TRUE for draws]\n")
  }
  cat("Formula:", deparse1(x$formula), "\n\n")

  # --- Variance components (lmer, fixed during sampling) ---
  cat("Random effects (variance components fixed at lmer estimates):\n")
  print(lme4::VarCorr(x$lmer), comp = "Std.Dev.", digits = digits)
  cat(sprintf("Number of obs: %d,  groups: %s, %d\n\n", n_obs, grp, n_grp))

  # --- Posterior means table ---
  cat("--- Posterior means (ICM exact, under fixed variance components) ---\n\n")

  # Flatten coef.mode to a data frame for aligned printing
  rows <- do.call(rbind, lapply(re_names, function(k) {
    nms <- names(x$coef.mode[[k]])
    data.frame(
      re  = k,
      par = nms,
      mode = unname(x$coef.mode[[k]]),
      stringsAsFactors = FALSE
    )
  }))

  w_re  <- max(nchar(rows$re),  nchar("RE component"))
  w_par <- max(nchar(rows$par), nchar("parameter"))

  if (!simulated) {
    # ICM-only: single value column
    cat(sprintf("  %-*s  %-*s  %12s\n",
                w_re, "RE component", w_par, "parameter", "coef.mode"))
    cat(sprintf("  %s  %s  %s\n",
                strrep("-", w_re), strrep("-", w_par), strrep("-", 12L)))
    for (i in seq_len(nrow(rows))) {
      cat(sprintf("  %-*s  %-*s  %12.*f\n",
                  w_re, rows$re[i], w_par, rows$par[i],
                  digits, rows$mode[i]))
    }
    cat("\n")

  } else {
    # Simulation: coef.mode + coef.means + draws SD side-by-side
    rows$means <- unlist(lapply(re_names, function(k) unname(x$coef.means[[k]])))
    rows$sd    <- unlist(lapply(re_names, function(k) {
      apply(x$fixef_draws[[k]], 2L, sd)
    }))

    cat(sprintf("  %-*s  %-*s  %12s  %12s  %10s\n",
                w_re, "RE component", w_par, "parameter",
                "coef.mode", "coef.means", "draws SD"))
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

  invisible(x)
}

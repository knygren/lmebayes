#' Bayesian generalized linear mixed model fit (draft)
#'
#' Draft entry point for \pkg{lmebayes} GLMM models with a \code{glmer}-like
#' interface, analogous to \code{\link{lmerb}} for Gaussian responses and to
#' \code{\link{glmb}} for fixed-effects GLMs.
#'
#' Currently a copy of \code{\link{lmerb}} with an additional \code{family}
#' argument. When \code{family = gaussian()}, behaviour matches
#' \code{\link{lmerb}} except that the embedded reference fit is from
#' \code{\link[lme4]{glmer}} rather than \code{\link[lme4]{lmer}}. Non-Gaussian
#' families use \code{block_rNormalGLM} for Block~1 Gibbs updates.
#'
#' @inheritParams lmerb
#' @param family A \code{\link[stats]{family}} object describing the response
#'   distribution and link. Defaults to \code{gaussian()}.
#' @param dispersion_ranef Observation-level measurement dispersion, treated
#'   as known during sampling.  Required positive scalar for families with a
#'   dispersion parameter (e.g. \code{gaussian()}); must be \code{NULL}
#'   (default) for \code{poisson()} and \code{binomial()}.  Typically
#'   \code{Prior_Setup_lmebayes(...)$dispersion_ranef}.
#' @param tv_tol Total variation tolerance per stored draw, in (0, 1)
#'   (default \code{0.01}).  For \code{family = gaussian()} the joint
#'   posterior is exactly multivariate normal and the number of inner Gibbs
#'   sweeps per stored draw is calibrated exactly as in \code{\link{lmerb}}
#'   (Nygren 2020, Theorem 3).  For non-Gaussian families the same
#'   calibration is applied to the \emph{local-Gaussian approximation of the
#'   posterior at its mode}: per-observation likelihood precisions are
#'   evaluated at the ICM posterior mode
#'   (\code{\link[glmbayesCore]{two_block_mode_weights}}) and fed to
#'   \code{\link[glmbayesCore]{two_block_rate_v2}}.  The derived sweep count is
#'   then the \emph{minimum} number of iterations required to converge to
#'   that hypothetical multivariate normal approximation -- a lower bound
#'   for the true (non-normal) posterior, not a guarantee.
#' @param m_convergence Optional integer override for the number of inner
#'   Gibbs sweeps per stored draw.  When \code{NULL} (default) the
#'   \code{tv_tol}-derived value is used.  A supplied value acts as a
#'   requested sweep count but is never allowed below the derived minimum:
#'   \code{max(m_convergence, m_min)} is used, with a warning if the value
#'   had to be raised.  Typical use is to pick a \emph{larger} number for
#'   non-Gaussian families (e.g. double the derived lower bound).
#' @param control Optional \code{\link[lme4]{glmerControl}} settings passed to
#'   the reference \code{\link[lme4]{glmer}} fit. Defaults to \code{NULL}
#'   (lme4 defaults). When \code{family = gaussian()}, lme4's \code{glmer}
#'   shortcut to \code{lmer} does not accept an explicit \code{glmerControl};
#'   leave \code{control = NULL} or pass \code{\link[lme4]{lmerControl}}.
#' @return Object of class \code{"glmerb"}: same structure as \code{"lmerb"},
#'   with additional \code{family} and \code{glmer} (reference
#'   \code{\link[lme4]{glmer}} fit) components instead of \code{lmer}.
#' @seealso \code{\link{lmerb}}, \code{\link[glmbayesCore]{glmerb_posterior_mode}},
#'   \code{\link{glmb}}
#' @examplesIf requireNamespace("bayesrules", quietly = TRUE)
#' @example inst/examples/Ex_glmerb.R
#' @export
glmerb <- function(
    formula,
    data = NULL,
    family = gaussian(),
    pfamily_list,
    dispersion_ranef = NULL,
    n = 1000L,
    tv_tol = 0.01,
    m_convergence = NULL,
    simulate = TRUE,
    REML = TRUE,
    control = NULL,
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
  if (missing(family) || is.null(family)) {
    family <- gaussian()
  }
  if (!inherits(family, "family")) {
    stop("'family' must be a family object.", call. = FALSE)
  }
  if (missing(pfamily_list) || is.null(pfamily_list)) {
    stop(
      "'pfamily_list' is required. Build it with ",
      "pfamily_list(Prior_Setup_lmebayes(...)) and pass the result to glmerb().",
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
    family = family,
    fit_mer = FALSE
  )

  design <- do.call(model_setup, c(setup_args, list(...)))
  if (!inherits(design, "model_setup")) {
    stop("model_setup() must return a model_setup object.", call. = FALSE)
  }

  prior <- .lmebayes_priors_from_pfamily_list(
    pfamily_list     = pfamily_list,
    dispersion_ranef = dispersion_ranef,
    design           = design,
    family           = family,
    fn_name          = "glmerb"
  )

  glmer_args <- c(
    list(
      formula = formula,
      data = data,
      family = family,
      verbose = verbose
    ),
    if (!is.null(control)) list(control = control),
    .lmebayes_mer_optional_args(
      start = start,
      subset = subset,
      weights = weights,
      na.action = na.action,
      offset = offset,
      contrasts = contrasts
    ),
    list(...)
  )
  glmer_fit <- do.call(lme4::glmer, glmer_args)

  if (is.null(fixef)) {
    fixef <- lapply(prior$prior_list, `[[`, "mu_fixef")
    names(fixef) <- design$re_coef_names
  }

  fixef_lmer <- fixef
  pm <- glmbayesCore::glmerb_posterior_mode(design, family, prior)
  fixef_start <- pm$fixef

  hdr <- sprintf("  %-18s  %-30s  %12s  %12s",
                 "RE component", "parameter", "lmer (start)", "post mode (ICM)")
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat("--- glmerb: Block 2 fixed effects ---\n")
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

  block1_prior <- .lmebayes_block1_prior_list(prior)

  if (!isTRUE(simulate)) {
    return(structure(
      list(
        call        = cl,
        formula     = formula,
        family      = family,
        glmer       = glmer_fit,
        prior       = prior,
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
      class = c("glmerb", "list")
    ))
  }

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

  # TV-calibrated number of inner Gibbs sweeps per stored draw.  Exact for
  # gaussian(), where the joint posterior is multivariate normal and the
  # Theorem 3 bound (Nygren 2020) applies.  For non-Gaussian families the
  # same machinery runs on the local-Gaussian approximation at the ICM
  # posterior mode (per-observation IRLS/Fisher weights at pm$b_mean), so the
  # derived m_min is the minimum number of sweeps required to converge to
  # that hypothetical multivariate normal approximation -- a lower bound for
  # the true posterior.  Chains start at the joint posterior mode (= the mean
  # of the approximating normal), so D0 = 0; + 1L covers the half-step lag of
  # the stored b draw.  A user-supplied m_convergence is floored at m_min.
  # ING components enter through the conservative disp_lower plug-in, making
  # lambda* an upper bound over the truncated tau^2 support.
  is_gaussian <- identical(family$family, "gaussian")
  if (is_gaussian) {
    rate <- glmbayesCore::two_block_rate_v2(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = block1_prior,
      pfamily_list      = prior$pfamily_list,
      family            = gaussian(),
      group_levels      = group_levels
    )
  } else {
    mode_w <- glmbayesCore::two_block_mode_weights(
      x            = design$Z,
      block        = design$groups,
      b_mode       = pm$b_mean,
      family       = family,
      group_levels = group_levels
    )
    rate <- glmbayesCore::two_block_rate_v2(
      x                 = design$Z,
      block             = design$groups,
      x_hyper           = design$X_hyper,
      prior_list_block1 = block1_prior,
      pfamily_list      = prior$pfamily_list,
      weights           = mode_w$weights,
      family            = family,
      group_levels      = group_levels
    )
  }
  m_min <- glmbayesCore::two_block_l_for_tv(
    rate, tv_tol, method = "theorem3"
  ) + 1L
  if (is.null(m_convergence)) {
    m_convergence <- m_min
  } else if (m_convergence < m_min) {
    warning(
      "glmerb: m_convergence = ", m_convergence, " is below the derived ",
      "minimum m_min = ", m_min, " for tv_tol = ", tv_tol,
      "; using m_min instead.",
      call. = FALSE
    )
    m_convergence <- m_min
  }
  calib_label <- if (is_gaussian) {
    "exact (Gaussian posterior)"
  } else {
    sprintf("approximate (local-Gaussian at mode, %s)", family$family)
  }
  if (prior$any_ing) {
    calib_label <- paste0(calib_label, "; conservative: ING tau^2_k = disp_lower")
  }
  cat(sprintf(
    "--- glmerb: convergence calibration [%s]:\n    lambda* = %.4f, tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n",
    calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
  ))
  method_label <- if (is_gaussian) "exact" else "local_gaussian_mode"
  if (prior$any_ing) {
    method_label <- paste0(method_label, "+disp_lower_bound")
  }
  convergence_info <- list(
    method        = method_label,
    tv_tol        = tv_tol,
    lambda_star   = rate$lambda_star,
    eigenvalues   = rate$eigenvalues,
    m_min         = m_min,
    m_convergence = m_convergence
  )

  # The v2 driver consumes the pfamily list directly: dNormal components get
  # the conjugate gamma_k draw at fixed tau^2_k (identical to the v1 path),
  # ING components make a joint (gamma_k, tau^2_k) draw via the
  # likelihood-subgradient envelope sampler, with the sampled tau^2_k fed
  # back into the Block 1 prior precision on the next inner step.
  sampler <- glmbayesCore::two_block_rNormal_reg_v2(
    n                 = n,
    y                 = design$y,
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = block1_prior,
    pfamily_list      = prior$pfamily_list,
    fixef_start       = fixef_start,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = design$group_name,
    family            = family,
    m_convergence     = m_convergence,
    seed              = seed,
    progbar           = TRUE
  )

  tau2_draws <- sampler$dispersion_fixef_draws

  structure(
    list(
      call         = cl,
      formula      = formula,
      family       = family,
      glmer        = glmer_fit,
      prior        = prior,
      model_setup  = design,
      coef.mode    = fixef_start,
      ranef.mode   = pm$b_mean,
      coef.means   = lapply(sampler$fixef_draws, colMeans),
      fixef_draws  = sampler$fixef_draws,
      coefficients = sampler$coefficients,
      tau2_draws   = tau2_draws,
      tau2.means   = colMeans(tau2_draws),
      mu_all       = sampler$mu_all_last,
      convergence  = convergence_info
    ),
    class = c("glmerb", "list")
  )
}

#' Print method for glmerb objects (draft)
#'
#' @param x Object of class \code{"glmerb"}.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.glmerb <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  re_names <- x$model_setup$re_coef_names
  grp      <- x$model_setup$group_name
  n_obs    <- length(x$model_setup$y)
  n_grp    <- nlevels(x$model_setup$groups)
  simulated <- !is.null(x$coefficients)
  fam      <- if (!is.null(x$family)) x$family$family else "gaussian"

  cat("Call:\n  ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  if (simulated) {
    n_draws <- nrow(x$fixef_draws[[re_names[1L]]])
    cat(sprintf(
      "Bayesian generalized linear mixed model  [%s; %d draws, two-block Gibbs]\n",
      fam, n_draws))
  } else {
    cat(sprintf(
      "Bayesian generalized linear mixed model  [%s; ICM only]\n", fam))
  }
  cat("Formula:", deparse1(x$formula), "\n\n")

  any_ing <- isTRUE(x$prior$any_ing)
  if (any_ing) {
    cat("Random effects (glmer reference; tau^2 sampled for ING components):\n")
  } else {
    cat("Random effects (variance components fixed at glmer estimates):\n")
  }
  print(lme4::VarCorr(x$glmer), comp = "Std.Dev.", digits = digits)
  cat(sprintf("Number of obs: %d,  groups: %s, %d\n\n", n_obs, grp, n_grp))
  if (any_ing && !is.null(x$tau2.means)) {
    cat("Posterior mean tau^2_k: ",
        paste(sprintf("%s = %.4g", names(x$tau2.means), x$tau2.means),
              collapse = ", "),
        "\n\n", sep = "")
  }

  cat("--- Posterior means (ICM exact, under fixed variance components) ---\n\n")

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

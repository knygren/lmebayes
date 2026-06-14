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
#' @param gap_tol Tolerated mode--mean gap in units of posterior standard
#'   deviations (default \code{0.0196}).  Applies only to non-Gaussian
#'   families.  The number of pilot chains is derived from this tolerance as
#'   \code{n_pilot = ceiling((qnorm(0.975) / gap_tol)^2)}, which ensures
#'   that a gap larger than \code{gap_tol} posterior SDs is detected with
#'   95\% power (two-sided, \eqn{\alpha = 0.05}).  The default
#'   \code{gap_tol = 0.0196 = 1.96 / 100} gives \code{n_pilot = 10000}.
#'   The pilot stage runs \code{n_pilot} independent chains from the ICM
#'   mode, each returning one stored draw after \code{m_convergence_pilot}
#'   inner sweeps, and takes the column-means as \code{coef.pilot.mean} --
#'   the starting point for the main run.  Set \code{gap_tol = NULL} to
#'   skip the pilot entirely and start the main run from the ICM mode (as
#'   in the Gaussian case).  Ignored for \code{family = gaussian()}, where
#'   mode equals mean exactly.
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
#' @param mode_gap_max Maximum per-coordinate mode--mean gap (in posterior
#'   standard deviation units) that \code{m_convergence_pilot} is calibrated
#'   to cover (default \code{1.0}).  Applies only to non-Gaussian families
#'   when \code{gap_tol} is not \code{NULL} and \code{m_convergence_pilot}
#'   is \code{NULL}.  The pilot chains start at the ICM mode, which is at
#'   Mahalanobis distance \eqn{D_{\max} = \sqrt{p}\,\times\,\texttt{mode\_gap\_max}}
#'   from the posterior mean (assuming \code{mode_gap_max} SDs per coordinate
#'   across \eqn{p} fixed-effect dimensions).  The number of pilot sweeps is
#'   the smallest \eqn{l} satisfying
#'   \eqn{\mathrm{erf}_1(0.5\,\lambda^{*l}\,D_{\max}/\sqrt{2}) \le
#'   \texttt{tv\_tol}} (Nygren 2020, Theorem 3 mean-shift term), floored at
#'   \code{m_min}.  Set \code{mode_gap_max = NULL} to fall back to
#'   \code{m_convergence} (pre-v0.2 behaviour).  Ignored for
#'   \code{family = gaussian()}.
#' @param m_convergence_pilot Optional integer override for the number of
#'   inner Gibbs sweeps used in each independent pilot chain.  Applies only
#'   when \code{gap_tol} is not \code{NULL} and \code{family} is
#'   non-Gaussian.  When \code{NULL} (default), the sweep count is derived
#'   automatically from \code{mode_gap_max}, \code{tv_tol}, and
#'   \code{rate$lambda_star} via Theorem 3 (see \code{mode_gap_max}).
#'   A supplied value is used as-is (no floor is applied beyond what the
#'   user intends) and overrides the \code{mode_gap_max} derivation.
#' @param control Optional \code{\link[lme4]{glmerControl}} settings passed to
#'   the reference \code{\link[lme4]{glmer}} fit. Defaults to \code{NULL}
#'   (lme4 defaults). When \code{family = gaussian()}, lme4's \code{glmer}
#'   shortcut to \code{lmer} does not accept an explicit \code{glmerControl};
#'   leave \code{control = NULL} or pass \code{\link[lme4]{lmerControl}}.
#' @return Object of class \code{"glmerb"}: same structure as \code{"lmerb"},
#'   with additional \code{family}, \code{glmer} (reference
#'   \code{\link[lme4]{glmer}} fit), \code{coef.pilot.mean} (estimated
#'   posterior mean from the pilot run, used as the main-run starting point;
#'   \code{NULL} for Gaussian or when \code{n_pilot = NULL}), and
#'   \code{pilot_mode_test} (multivariate Wald test of pilot mean against
#'   ICM mode: \code{Q}, \code{df}, \code{p_value},
#'   \code{n_pilot}; \code{NULL} when no pilot), \code{gap_tol}
#'   (the tolerance used to derive \code{n_pilot}), and \code{mode_gap_max}
#'   (the per-coordinate gap tolerance used to derive
#'   \code{m_convergence_pilot}) components instead of \code{lmer}.
#' @seealso \code{\link{lmerb}}, \code{\link[glmbayesCore]{glmerb_posterior_mode}},
#'   \code{\link{glmb}}; \code{\link[utils]{demo}} for the full sampling workflow
#'   (\code{demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")}).
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
    gap_tol = 0.0196,
    mode_gap_max = 1.0,
    tv_tol = 0.01,
    m_convergence = NULL,
    m_convergence_pilot = NULL,
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
  if (!is.null(gap_tol)) {
    if (!is.numeric(gap_tol) || length(gap_tol) != 1L ||
        !is.finite(gap_tol) || gap_tol <= 0 || gap_tol >= 1) {
      stop("'gap_tol' must be NULL or a single value in (0, 1).", call. = FALSE)
    }
  }
  n_pilot <- if (!is.null(gap_tol)) {
    as.integer(ceiling((stats::qnorm(0.975) / gap_tol)^2))
  } else {
    NULL
  }
  if (!is.null(mode_gap_max)) {
    if (!is.numeric(mode_gap_max) || length(mode_gap_max) != 1L ||
        !is.finite(mode_gap_max) || mode_gap_max <= 0) {
      stop("'mode_gap_max' must be NULL or a single positive finite number.",
           call. = FALSE)
    }
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
  if (!is.null(m_convergence_pilot)) {
    if (!is.numeric(m_convergence_pilot) ||
        length(m_convergence_pilot) != 1L ||
        !is.finite(m_convergence_pilot) ||
        m_convergence_pilot < 1) {
      stop("'m_convergence_pilot' must be NULL or a single integer >= 1.",
           call. = FALSE)
    }
    m_convergence_pilot <- as.integer(m_convergence_pilot)
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

  if (!isTRUE(simulate)) {
    fixef_glmer <- fixef
    pm          <- glmbayesCore::glmerb_posterior_mode(design, family, prior)
    fixef_start <- pm$fixef
    hdr <- sprintf("  %-18s  %-30s  %12s  %12s",
                   "RE component", "parameter", "glmer (start)", "post mode (ICM)")
    sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
    cat("--- glmerb: Block 2 fixed effects ---\n")
    cat(hdr, "\n")
    cat(sep, "\n")
    for (k in design$re_coef_names) {
      nms_k   <- names(fixef_glmer[[k]])
      glmer_v <- fixef_glmer[[k]]
      pm_v    <- fixef_start[[k]]
      for (nm in nms_k) {
        cat(sprintf("  %-18s  %-30s  %12.4f  %12.4f\n",
                    k, nm, glmer_v[[nm]], pm_v[[nm]]))
      }
    }
    cat(sprintf("  (ICM converged: %s, %d iter, delta = %.2e)\n\n",
                pm$converged, pm$iterations, pm$delta))
    return(structure(
      list(
        call         = cl,
        formula      = formula,
        family       = family,
        glmer        = glmer_fit,
        prior        = prior,
        model_setup  = design,
        coef.mode    = fixef_start,
        ranef.mode   = pm$b_mean,
        coef.means   = NULL,
        fixef_draws  = NULL,
        coefficients = NULL,
        mu_all       = as.matrix(
          glmbayesCore::build_mu_all(design, fixef_start)$mu_all
        )
      ),
      class = c("glmerb", "list")
    ))
  }

  run_pilot <- !identical(family$family, "gaussian") && !is.null(n_pilot)

  # ICM mode, convergence calibration, pilot stage, and main stage are all
  # handled inside rglmerb.  glmerb passes design + prior so that rglmerb can
  # call glmerb_posterior_mode() internally.
  sampler <- rglmerb(
    n                   = n,
    design              = design,
    prior               = prior,
    family              = family,
    fixef_start         = NULL,           # computed internally from design + prior
    m_convergence       = m_convergence,  # NULL => derived from tv_tol
    n_pilot             = if (run_pilot) n_pilot else 0L,
    m_convergence_pilot = m_convergence_pilot,
    tv_tol              = tv_tol,
    mode_gap_max        = mode_gap_max,
    seed                = seed,
    collect_block1      = TRUE,
    verbose             = TRUE
  )

  convergence_info <- sampler$convergence_info
  pilot_mode_test  <- sampler$pilot_mode_test
  fixef_start      <- sampler$coef.mode
  fixef_main_start <- sampler$fixef_main_start
  tau2_draws       <- sampler$dispersion_fixef_draws
  iters_draws      <- sampler$iters_fixef_draws
  m_convergence    <- sampler$m_convergence_used

  structure(
    list(
      call              = cl,
      formula           = formula,
      family            = family,
      glmer             = glmer_fit,
      prior             = prior,
      model_setup       = design,
      coef.mode         = fixef_start,
      coef.pilot.mean   = if (run_pilot) fixef_main_start else NULL,
      ranef.mode        = sampler$ranef.mode,
      coef.means        = lapply(sampler$fixef_draws, colMeans),
      fixef_draws       = sampler$fixef_draws,
      coefficients      = sampler$coefficients,
      tau2_draws        = tau2_draws,
      tau2.means        = colMeans(tau2_draws),
      iters_draws       = iters_draws,
      iters.means       = colMeans(iters_draws) / m_convergence,
      mu_all            = sampler$mu_all_last,
      pilot_mode_test   = pilot_mode_test,
      gap_tol           = gap_tol,
      mode_gap_max      = mode_gap_max,
      convergence       = convergence_info
    ),
    class = c("glmerb", "list")
  )
}

#' Print posterior estimates by RE component for a glmerb / lmerb fit
#'
#' Displays a side-by-side table of the lmer/glmer MLE reference, posterior
#' mode or ICM mean (\code{coef.mode}), and posterior mean (\code{coef.means})
#' for every (RE-component, parameter) pair.  When \code{x} is a bare
#' \code{coef.means} list rather than a full fit object, only the posterior
#' mean column is shown.
#'
#' For \code{lmerb} objects the reference column is labelled \code{"lmer"} and
#' the \code{coef.mode} column is labelled \code{"ICM.mean"} (the Gaussian
#' posterior mean and mode coincide exactly).  For \code{glmerb} objects the
#' reference column is labelled \code{"glmer"} and the \code{coef.mode} column
#' is labelled \code{"post.mode"}.
#'
#' @param x A \code{glmerb} or \code{lmerb} object, or a bare \code{coef.means}
#'   list.
#' @param digits Number of decimal places for numeric columns.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print_coef_means <- function(x, digits = 4L, ...) {
  is_fit    <- inherits(x, c("glmerb", "lmerb"))
  is_lmerb  <- inherits(x, "lmerb")
  cm        <- if (is_fit) x$coef.means else x
  if (is.null(cm)) {
    cat("coef.means: NULL (simulation not yet run)\n")
    return(invisible(x))
  }

  rows <- do.call(rbind, lapply(names(cm), function(k) {
    v <- cm[[k]]
    data.frame(component = k, parameter = names(v),
               post_mean = unname(v), stringsAsFactors = FALSE)
  }))

  # Reference MLE column: glmer for glmerb, lmer for lmerb.
  mer_fit   <- if (is_fit) (if (is_lmerb) x$lmer else x$glmer) else NULL
  has_mer   <- !is.null(mer_fit)
  mer_label <- if (is_lmerb) "lmer" else "glmer"

  if (has_mer) {
    mer_v <- lme4::fixef(mer_fit)
    # Map (component, parameter) -> fixef name using the same convention as
    # fe_name_for() in Prior_Setup_lmebayes:
    #   (Intercept) component, col X  -> fixef["X"]
    #   component K, (Intercept) col  -> fixef["K"]
    #   component K, col X            -> fixef["X:K"] or fixef["K:X"]
    rows[[mer_label]] <- mapply(function(k, col) {
      nm <- if (k == "(Intercept)") {
        col
      } else if (col == "(Intercept)") {
        k
      } else {
        cand <- c(paste0(col, ":", k), paste0(k, ":", col))
        hit  <- cand[cand %in% names(mer_v)]
        if (length(hit)) hit[1L] else NA_character_
      }
      if (!is.na(nm) && nm %in% names(mer_v)) unname(mer_v[nm]) else NA_real_
    }, rows$component, rows$parameter)
  }

  # coef.mode column: "ICM.mean" for lmerb (exact posterior mean), "post.mode"
  # for glmerb (posterior mode from ICM optimisation).
  has_mode   <- is_fit && !is.null(x$coef.mode)
  mode_label <- if (is_lmerb) "ICM.mean" else "post.mode"
  if (has_mode) {
    rows[[mode_label]] <- unlist(lapply(names(x$coef.mode), function(k) {
      unname(x$coef.mode[[k]])
    }))
  }

  w_c <- max(nchar(rows$component), nchar("RE component"))
  w_p <- max(nchar(rows$parameter),  nchar("parameter"))
  w_v <- digits + 4L

  cols <- character(0)
  if (has_mer)  cols <- c(cols, mer_label)
  if (has_mode) cols <- c(cols, mode_label)
  cols <- c(cols, "post.mean")
  n_val <- length(cols)

  # Pre-format numeric values so the outer sprintf only ever sees %s.
  num_fmt <- sprintf("%%%d.%df", w_v, digits)   # e.g. "%8.4f"
  val_hdr <- paste(formatC(cols, width = w_v, flag = " "), collapse = "  ")
  val_sep <- paste(rep(strrep("-", w_v), n_val), collapse = "  ")

  cat(sprintf("  %-*s  %-*s  %s\n", w_c, "RE component", w_p, "parameter", val_hdr))
  cat(sprintf("  %-*s  %-*s  %s\n", w_c, strrep("-", w_c), w_p, strrep("-", w_p), val_sep))
  for (i in seq_len(nrow(rows))) {
    vals <- character(0L)
    if (has_mer)  vals <- c(vals, sprintf(num_fmt, rows[[mer_label]][i]))
    if (has_mode) vals <- c(vals, sprintf(num_fmt, rows[[mode_label]][i]))
    vals <- c(vals, sprintf(num_fmt, rows$post_mean[i]))
    cat(sprintf("  %-*s  %-*s  %s\n",
                w_c, rows$component[i],
                w_p, rows$parameter[i],
                paste(vals, collapse = "  ")))
  }
  invisible(x)
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

#' Bayesian generalized linear mixed model fit (draft)
#'
#' Draft entry point for \pkg{lmebayes} GLMM models with a \code{glmer}-like
#' interface, analogous to \code{\link{lmerb}} for Gaussian responses and to
#' \code{\link{glmb}} for fixed-effects GLMs.
#'
#' Currently a copy of \code{\link{lmerb}} with an additional \code{family}
#' argument. Only \code{family = gaussian()} is implemented; other families
#' stop with an informative error until GLMM Block~1 samplers are wired in.
#'
#' @inheritParams lmerb
#' @param family A \code{\link[stats]{family}} object describing the response
#'   distribution and link. Defaults to \code{gaussian()}. When
#'   \code{gaussian()}, behaviour matches \code{\link{lmerb}}.
#' @return Object of class \code{"glmerb"}: same structure as \code{"lmerb"},
#'   with an additional \code{family} component.
#' @seealso \code{\link{lmerb}}, \code{\link[glmbayesCore]{glmerb_posterior_mode}},
#'   \code{\link{glmb}}
#' @examples
#' \donttest{
#'   source(system.file("examples", "Ex_glmerb.R", package = "lmebayes"))
#' }
#' @export
glmerb <- function(
    formula,
    data = NULL,
    family = gaussian(),
    measurement_prior_list,
    n = 1000L,
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
  if (missing(family) || is.null(family)) {
    family <- gaussian()
  }
  if (!inherits(family, "family")) {
    stop("'family' must be a family object.", call. = FALSE)
  }
  if (!identical(family$family, "gaussian")) {
    stop(
      "glmerb() is a draft: only family = gaussian() is implemented. ",
      "Use lmerb() for the Gaussian mixed model.",
      call. = FALSE
    )
  }
  if (missing(measurement_prior_list) || is.null(measurement_prior_list)) {
    stop(
      "'measurement_prior_list' is required. ",
      "Build it with Prior_Setup_lmebayes() and pass the result to glmerb().",
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

  fixef_lmer <- fixef
  pm <- glmbayesCore::glmerb_posterior_mode(design, family, measurement_prior_list)
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

  Sigma_ranef <- measurement_prior_list$Sigma_ranef
  if (is.null(Sigma_ranef)) {
    stop("measurement_prior_list must contain 'Sigma_ranef'.", call. = FALSE)
  }
  P <- solve(Sigma_ranef)

  dispersion <- measurement_prior_list$dispersion_ranef
  if (is.null(dispersion)) {
    stop("measurement_prior_list must contain 'dispersion_ranef'.", call. = FALSE)
  }

  if (!isTRUE(simulate)) {
    return(structure(
      list(
        call        = cl,
        formula     = formula,
        family      = family,
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
      class = c("glmerb", "list")
    ))
  }

  m_convergence <- 10L

  re_names     <- design$re_coef_names
  group_levels <- levels(design$groups)

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

  sampler <- glmbayesCore::two_block_rNormal_reg(
    n                 = n,
    y                 = design$y,
    x                 = design$Z,
    block             = design$groups,
    x_hyper           = design$X_hyper,
    prior_list_block1 = list(P = P, dispersion = dispersion, ddef = FALSE),
    prior_list_block2 = block2_prior_list,
    fixef_start       = fixef_start,
    re_coef_names     = re_names,
    group_levels      = group_levels,
    group_name        = design$group_name,
    m_convergence     = m_convergence,
    seed              = seed,
    progbar           = TRUE
  )

  structure(
    list(
      call         = cl,
      formula      = formula,
      family       = family,
      lmer         = lmer_fit,
      prior        = measurement_prior_list,
      model_setup  = design,
      coef.mode    = fixef_start,
      ranef.mode   = pm$b_mean,
      coef.means   = lapply(sampler$fixef_draws, colMeans),
      fixef_draws  = sampler$fixef_draws,
      coefficients = sampler$coefficients,
      mu_all       = sampler$mu_all_last
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

  cat("Random effects (variance components fixed at lmer estimates):\n")
  print(lme4::VarCorr(x$lmer), comp = "Std.Dev.", digits = digits)
  cat(sprintf("Number of obs: %d,  groups: %s, %d\n\n", n_obs, grp, n_grp))

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

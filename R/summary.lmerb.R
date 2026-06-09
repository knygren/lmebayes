#' Summarize Bayesian linear mixed model fits
#'
#' @description
#' Methods for \code{\link{lmerb}} fits.  \code{summary.lmerb} builds
#' Block~2 (level-2 fixed effect) tables per random-effects component,
#' following the layout of \code{\link[glmbayes]{summary.glmb}} and the
#' multi-response structure of \code{\link[glmbayes]{summary.mlmb}}.
#'
#' @param object An object of class \code{"lmerb"} or \code{"glmerb"}.
#' @param groups Optional character vector of grouping levels for which to
#'   include a per-group Block~1 (random effects) detail table.  When
#'   \code{NULL} (default), only an aggregate \code{ranef_overview} is
#'   returned.
#' @param digits Number of significant digits for printing.
#' @param \ldots Ignored.
#' @return \code{summary.lmerb} returns an object of class
#'   \code{"summary.lmerb"}, a list with components \code{call},
#'   \code{formula}, \code{n}, \code{simulated}, \code{varcor},
#'   \code{fixef_overview}, \code{fixef} (per-RE-component tables),
#'   \code{ranef_overview}, and optionally \code{ranef_groups}.
#' @seealso \code{\link{lmerb}}, \code{\link{glmerb}}, \code{\link{print.lmerb}},
#'   \code{\link[glmbayes]{summary.glmb}}, \code{\link[glmbayes]{summary.mlmb}}
#' @export
#' @method summary lmerb
summary.lmerb <- function(object, groups = NULL, digits = max(3L, getOption("digits") - 3L), ...) {

  if (!inherits(object, c("lmerb", "glmerb"))) {
    stop("'object' must be an lmerb or glmerb fit.", call. = FALSE)
  }

  re_names  <- object$model_setup$re_coef_names
  simulated <- !is.null(object$coefficients)
  n_draws   <- if (simulated) nrow(object$fixef_draws[[re_names[1L]]]) else NULL

  fixef_parts <- stats::setNames(
    lapply(re_names, function(k) {
      .lmerb_fixef_component_summary(object, k, n_draws = n_draws, simulated = simulated)
    }),
    re_names
  )

  res <- list(
    call          = object$call,
    formula       = object$formula,
    n             = n_draws,
    simulated     = simulated,
    lmer          = object$lmer,
    varcor        = lme4::VarCorr(object$lmer),
    dispersion    = object$prior$dispersion_ranef,
    n_obs         = length(object$model_setup$y),
    n_groups      = nlevels(object$model_setup$groups),
    group_name    = object$model_setup$group_name,
    fixef_overview = .lmerb_fixef_overview(object, simulated = simulated),
    fixef         = fixef_parts,
    ranef_overview = .lmerb_ranef_overview(object, simulated = simulated)
  )

  if (!is.null(groups) && length(groups) > 0L) {
    res$ranef_groups <- .lmerb_ranef_groups_detail(object, groups, simulated = simulated)
  }

  class(res) <- "summary.lmerb"
  res
}

#' @rdname summary.lmerb
#' @export
#' @method summary glmerb
summary.glmerb <- summary.lmerb

#' @rdname summary.lmerb
#' @param x An object of class \code{"summary.lmerb"}.
#' @export
#' @method print summary.lmerb
print.summary.lmerb <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("Call:\n  ")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"))
  cat("\n\n")

  if (isTRUE(x$simulated)) {
    cat(sprintf("Bayesian linear mixed model fit  [%d draws, two-block Gibbs]\n", x$n))
  } else {
    cat("Bayesian linear mixed model fit  [ICM only; simulation not run]\n")
  }
  cat("Formula:", deparse1(x$formula), "\n\n")

  cat("Random effects (variance components fixed at lmer estimates):\n")
  print(x$varcor, comp = "Std.Dev.", digits = digits)
  cat(sprintf(
    "Number of obs: %d,  groups: %s, %d\n\n",
    x$n_obs, x$group_name, x$n_groups
  ))

  # --- Block 2 overview ---
  cat("=== Block 2: Level-2 fixed effects (hyperparameters) ===\n\n")
  cat("Overview:\n")
  if (!is.null(x$fixef_overview) && nrow(x$fixef_overview) > 0L) {
    stats::printCoefmat(x$fixef_overview, digits = digits, quote = FALSE)
  } else {
    cat("  (no fixed-effect hyperparameters)\n")
  }
  cat("\n")

  # --- Per RE component ---
  for (k in names(x$fixef)) {
    part <- x$fixef[[k]]
    cat("--- RE component:", k, "---\n\n")

    cat("Prior and lmer reference:\n\n")
    if (!is.null(part$coefficients1)) {
      stats::printCoefmat(part$coefficients1, digits = digits, quote = FALSE)
    }
    cat("\n")

    if (isTRUE(x$simulated)) {
      cat("Bayesian estimates based on", x$n, "draws:\n\n")
      if (!is.null(part$coefficients)) {
        stats::printCoefmat(part$coefficients, digits = digits, quote = FALSE)
      }
      cat("\nDistribution percentiles:\n\n")
      if (!is.null(part$Percentiles)) {
        stats::printCoefmat(part$Percentiles, digits = digits, quote = FALSE)
      }
    } else {
      cat("Posterior mode (= mean, ICM exact):\n\n")
      if (!is.null(part$coefficients)) {
        stats::printCoefmat(
          part$coefficients[, "Post.Mode", drop = FALSE],
          digits = digits,
          quote = FALSE
        )
      }
      cat("\n  (Run with simulate = TRUE for MCMC means, SDs, and percentiles.)\n")
    }
    cat("\n")
  }

  # --- Block 1 overview ---
  cat("=== Block 1: Random effects (group-level) ===\n\n")
  cat("Summary of posterior mode (ranef.mode) across groups:\n\n")
  if (!is.null(x$ranef_overview) && nrow(x$ranef_overview) > 0L) {
    stats::printCoefmat(x$ranef_overview, digits = digits, quote = FALSE)
  }
  cat("\n")

  if (!is.null(x$ranef_groups)) {
    cat("Per-group detail (requested levels):\n\n")
    print(round(x$ranef_groups, digits))
    cat("\n")
  } else {
    cat(
      "Per-group random effects: inspect fit$ranef.mode or fit$coefficients,\n",
      "  or call summary(fit, groups = <level ids>) for selected groups.\n\n",
      sep = ""
    )
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
.lmerb_lmer_fixef_lookup <- function(lmer_fit, re_name, par_name) {
  fe <- lme4::fixef(lmer_fit)
  fe_names <- names(fe)

  candidates <- character(0)
  if (par_name == "(Intercept)" && re_name == "(Intercept)") {
    candidates <- c("(Intercept)")
  } else if (par_name == "(Intercept)" && re_name != "(Intercept)") {
    candidates <- c(re_name)
  } else {
    candidates <- c(
      par_name,
      paste0(par_name, ":", re_name),
      paste0(re_name, ":", par_name)
    )
  }

  hit <- candidates[candidates %in% fe_names]
  if (length(hit) == 0L) {
    return(list(estimate = NA_real_, se = NA_real_))
  }

  nm <- hit[1L]
  sm <- tryCatch(summary(lmer_fit), error = function(e) NULL)
  est <- unname(fe[[nm]])
  se  <- NA_real_
  if (!is.null(sm) && nm %in% rownames(sm$coefficients)) {
    se <- sm$coefficients[nm, "Std. Error"]
  }
  list(estimate = est, se = se)
}

#' @keywords internal
.lmerb_fixef_component_summary <- function(object, k, n_draws, simulated) {

  pl_k   <- object$prior$prior_list[[k]]
  par    <- names(object$coef.mode[[k]])
  q_k    <- length(par)

  prior_mean <- unname(pl_k$mu_fixef)
  prior_sd   <- sqrt(diag(pl_k$Sigma_fixef))

  lmer_ref <- lapply(par, function(nm) {
    .lmerb_lmer_fixef_lookup(object$lmer, k, nm)
  })
  lmer_est <- vapply(lmer_ref, `[[`, numeric(1), "estimate")
  lmer_se  <- vapply(lmer_ref, `[[`, numeric(1), "se")

  Tab1 <- cbind(
    "Prior Mean" = prior_mean,
    "Prior.sd"   = prior_sd,
    "lmer"       = lmer_est,
    "lmer.se"    = lmer_se
  )
  rownames(Tab1) <- par

  post_mode <- unname(object$coef.mode[[k]])

  if (!simulated) {
    TAB <- cbind("Post.Mode" = post_mode)
    rownames(TAB) <- par
    return(list(
      coefficients1 = Tab1,
      coefficients  = TAB,
      Percentiles   = NULL
    ))
  }

  draws <- object$fixef_draws[[k]]
  post_mean <- unname(object$coef.means[[k]])
  post_sd   <- apply(draws, 2L, stats::sd)
  mc_err    <- post_sd / sqrt(n_draws)

  pval2 <- vapply(seq_len(q_k), function(j) {
    p1 <- mean(draws[, j] < prior_mean[j])
    min(p1, 1 - p1)
  }, numeric(1))

  percentiles <- t(apply(draws, 2L, stats::quantile,
    probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
  ))

  TAB <- cbind(
    "Post.Mode"      = post_mode,
    "Post.Mean"      = post_mean,
    "Post.Sd"        = post_sd,
    "MC Error"       = mc_err,
    "Pr(Prior_tail)" = pval2
  )
  rownames(TAB) <- par

  TAB2 <- cbind(
    "1.0%"  = percentiles[, 1],
    "2.5%"  = percentiles[, 2],
    "5.0%"  = percentiles[, 3],
    "Median" = percentiles[, 4],
    "95.0%" = percentiles[, 5],
    "97.5%" = percentiles[, 6],
    "99.0%" = percentiles[, 7]
  )
  rownames(TAB2) <- par

  list(
    coefficients1 = Tab1,
    coefficients  = TAB,
    Percentiles   = TAB2
  )
}

#' @keywords internal
.lmerb_fixef_overview <- function(object, simulated) {

  re_names <- object$model_setup$re_coef_names
  rows <- do.call(rbind, lapply(re_names, function(k) {
    nms <- names(object$coef.mode[[k]])
    data.frame(
      parameter = paste0(k, "::", nms),
      coef.mode = unname(object$coef.mode[[k]]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(rows) <- rows$parameter
  out <- rows[, "coef.mode", drop = FALSE]
  colnames(out) <- "coef.mode"

  if (simulated) {
    means <- unlist(lapply(re_names, function(k) unname(object$coef.means[[k]])))
    sds   <- unlist(lapply(re_names, function(k) {
      apply(object$fixef_draws[[k]], 2L, stats::sd)
    }))
    n     <- nrow(object$fixef_draws[[re_names[1L]]])
    out   <- cbind(
      out,
      coef.means = means,
      Post.Sd    = sds,
      MC.Error   = sds / sqrt(n)
    )
  }

  out
}

#' @keywords internal
.lmerb_ranef_overview <- function(object, simulated) {

  re_names <- object$model_setup$re_coef_names
  b_mode   <- object$ranef.mode

  overview <- t(vapply(re_names, function(k) {
    v <- b_mode[, k]
    c(
      Mean   = mean(v),
      SD     = stats::sd(v),
      Min    = min(v),
      Q1     = unname(stats::quantile(v, 0.25)),
      Median = median(v),
      Q3     = unname(stats::quantile(v, 0.75)),
      Max    = max(v)
    )
  }, numeric(7)))

  if (simulated) {
    grp_col <- object$model_setup$group_name
    mcmc_means <- vapply(re_names, function(k) {
      mean(tapply(
        object$coefficients[[k]],
        object$coefficients[[grp_col]],
        mean
      ))
    }, numeric(1))
    overview <- cbind(overview, MCMC.mean = mcmc_means)
  }

  overview
}

#' @keywords internal
.lmerb_ranef_groups_detail <- function(object, groups, simulated) {

  re_names <- object$model_setup$re_coef_names
  grp_col  <- object$model_setup$group_name
  groups   <- as.character(groups)

  rows <- lapply(groups, function(lev) {
    mode_vals <- object$ranef.mode[lev, re_names, drop = TRUE]
    out <- data.frame(
      group = lev,
      t(mode_vals),
      check.names = FALSE
    )
    names(out) <- c("group", re_names)

    if (simulated) {
      idx <- object$coefficients[[grp_col]] == lev
      for (k in re_names) {
        out[[paste0(k, ".mean")]] <- mean(object$coefficients[idx, k])
        out[[paste0(k, ".sd")]]   <- stats::sd(object$coefficients[idx, k])
      }
    }
    out
  })

  do.call(rbind, rows)
}

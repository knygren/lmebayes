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
#'   \code{fixef_prior_overview} (stacked prior and \code{glmer}/\code{lmer}
#'   reference across RE components),
#'   \code{fixef_overview} (Block~2 hyperparameters with posterior summaries and
#'   \code{Pr(Prior_tail)}),
#'   \code{fixef_percentiles_overview} (stacked distribution percentiles across
#'   RE components, when simulated),
#'   \code{fixef} (per-RE-component tables; not printed, available on the
#'   returned object),
#'   \code{ranef_overview}, \code{ranef.iters.mean} (Block~1 envelope candidates
#'   per inner sweep, averaged over groups; printed separately from the overview
#'   table), \code{any_non_normal}, \code{tau2} (per-component Block~2 dispersion
#'   prior type, plug-in value, posterior mean / SD / quantiles from
#'   \code{fixef.dispersion} for sampled ING components, and \code{Cand/draw},
#'   the average number of Block~2 candidates per accepted draw from
#'   \code{fixef.iters.mean} -- 1 for \code{dNormal}, roughly the reciprocal
#'   acceptance rate of the envelope sampler for ING), and optionally
#'   \code{ranef_groups}.
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
  n_draws   <- if (simulated) nrow(object$fixef[[re_names[1L]]]) else NULL
  mer_fit   <- .lmerb_reference_fit(object)
  mer_label <- if (inherits(object, "glmerb")) "glmer" else "lmer"

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
    mer_label     = mer_label,
    mer           = mer_fit,
    varcor        = lme4::VarCorr(mer_fit),
    dispersion    = object$prior$dispersion_ranef,
    n_obs         = length(object$model_setup$y),
    n_groups      = nlevels(object$model_setup$groups),
    group_name    = object$model_setup$group_name,
    fixef_prior_overview = .lmerb_fixef_prior_overview(fixef_parts),
    fixef_overview = .lmerb_fixef_overview(object, simulated = simulated),
    fixef_percentiles_overview = .lmerb_fixef_percentiles_overview(fixef_parts),
    fixef         = fixef_parts,
    ranef_overview = .lmerb_ranef_overview(object, simulated = simulated),
    any_non_normal = isTRUE(object$any_non_normal) ||
      isTRUE(object$prior$any_non_normal),
    tau2          = .lmerb_tau2_summary(object, simulated = simulated),
    ranef.iters.mean = if (simulated) object$ranef.iters.mean else NULL
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

  mer_label <- if (!is.null(x$mer_label)) x$mer_label else "lmer"
  if (isTRUE(x$any_non_normal)) {
    cat(sprintf(
      "Random effects (%s reference; tau^2 sampled for non-dNormal components):\n",
      mer_label
    ))
  } else {
    cat(sprintf(
      "Random effects (variance components fixed at %s estimates):\n",
      mer_label
    ))
  }
  print(x$varcor, comp = "Std.Dev.", digits = digits)
  cat(sprintf(
    "Number of obs: %d,  groups: %s, %d\n\n",
    x$n_obs, x$group_name, x$n_groups
  ))

  if (!is.null(x$tau2)) {
    cat("Block 2 dispersion (RE variance tau^2_k):\n\n")
    t2 <- x$tau2
    num_cols <- vapply(t2, is.numeric, logical(1L))
    t2[num_cols] <- lapply(t2[num_cols], round, digits = digits)
    print(t2)
    cat("\n")
  }

  # --- Block 2 overview ---
  cat("=== Block 2: Level-2 fixed effects (hyperparameters) ===\n\n")
  cat(sprintf("Prior and %s reference:\n\n", mer_label))
  if (!is.null(x$fixef_prior_overview) && nrow(x$fixef_prior_overview) > 0L) {
    stats::printCoefmat(x$fixef_prior_overview, digits = digits, quote = FALSE)
  } else {
    cat("  (no fixed-effect hyperparameters)\n")
  }
  cat("\nOverview:\n")
  if (!is.null(x$fixef_overview) && nrow(x$fixef_overview) > 0L) {
    stats::printCoefmat(x$fixef_overview, digits = digits, quote = FALSE)
  } else {
    cat("  (no fixed-effect hyperparameters)\n")
  }
  if (isTRUE(x$simulated)) {
    cat("\nDistribution percentiles:\n\n")
    if (!is.null(x$fixef_percentiles_overview) &&
        nrow(x$fixef_percentiles_overview) > 0L) {
      stats::printCoefmat(x$fixef_percentiles_overview, digits = digits, quote = FALSE)
    }
  } else {
    cat("\n  (Run with simulate = TRUE for MCMC means, SDs, and percentiles.)\n")
  }
  cat("\n")

  # --- Block 1 overview ---
  cat("=== Block 1: Random effects (group-level) ===\n\n")
  cat("Summary of posterior mode (ranef.mode) across groups:\n\n")
  if (!is.null(x$ranef_overview) && nrow(x$ranef_overview) > 0L) {
    stats::printCoefmat(x$ranef_overview, digits = digits, quote = FALSE)
  }
  if (isTRUE(x$simulated) && !is.null(x$ranef.iters.mean)) {
    cat(
      "\nMean Block 1 likelihood subgradient candidates per stored draw:",
      formatC(x$ranef.iters.mean, digits = digits, format = "f"),
      "\n  (averaged over groups; same for all RE components in a sweep)\n\n"
    )
  } else {
    cat("\n")
  }

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
  par    <- names(object$fixef.mode[[k]])
  q_k    <- length(par)

  prior_mean <- unname(pl_k$mu_fixef)
  prior_sd   <- sqrt(diag(pl_k$Sigma_fixef))

  mer_ref <- lapply(par, function(nm) {
    .lmerb_lmer_fixef_lookup(.lmerb_reference_fit(object), k, nm)
  })
  mer_est <- vapply(mer_ref, `[[`, numeric(1), "estimate")
  mer_se  <- vapply(mer_ref, `[[`, numeric(1), "se")
  mer_label <- if (inherits(object, "glmerb")) "glmer" else "lmer"

  Tab1 <- cbind(
    "Prior Mean" = prior_mean,
    "Prior.sd"   = prior_sd,
    mer_est,
    mer_se
  )
  colnames(Tab1) <- c("Prior Mean", "Prior.sd", mer_label, paste0(mer_label, ".se"))
  rownames(Tab1) <- par

  post_mode <- unname(object$fixef.mode[[k]])

  if (!simulated) {
    TAB <- cbind("Post.Mode" = post_mode)
    rownames(TAB) <- par
    return(list(
      coefficients1 = Tab1,
      coefficients  = TAB,
      Percentiles   = NULL
    ))
  }

  draws <- object$fixef[[k]]
  post_mean <- unname(object$fixef.means[[k]])
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
.lmerb_fixef_prior_overview <- function(fixef_parts) {

  if (length(fixef_parts) == 0L) {
    return(NULL)
  }

  rows_list <- lapply(names(fixef_parts), function(k) {
    tab <- fixef_parts[[k]]$coefficients1
    if (is.null(tab) || nrow(tab) == 0L) {
      return(NULL)
    }
    rownames(tab) <- paste0(k, "::", rownames(tab))
    tab
  })
  rows_list <- rows_list[!vapply(rows_list, is.null, logical(1L))]
  if (length(rows_list) == 0L) {
    return(NULL)
  }

  do.call(rbind, rows_list)
}

#' @keywords internal
.lmerb_fixef_percentiles_overview <- function(fixef_parts) {

  if (length(fixef_parts) == 0L) {
    return(NULL)
  }

  rows_list <- lapply(names(fixef_parts), function(k) {
    tab <- fixef_parts[[k]]$Percentiles
    if (is.null(tab) || nrow(tab) == 0L) {
      return(NULL)
    }
    rownames(tab) <- paste0(k, "::", rownames(tab))
    tab
  })
  rows_list <- rows_list[!vapply(rows_list, is.null, logical(1L))]
  if (length(rows_list) == 0L) {
    return(NULL)
  }

  do.call(rbind, rows_list)
}

#' @keywords internal
.lmerb_fixef_overview <- function(object, simulated) {

  re_names <- object$model_setup$re_coef_names

  rows_list <- lapply(re_names, function(k) {
    par <- names(object$fixef.mode[[k]])
    n_p <- length(par)

    out <- data.frame(
      fixef.mode = unname(object$fixef.mode[[k]]),
      stringsAsFactors = FALSE
    )

    if (simulated) {
      draws      <- object$fixef[[k]]
      post_mean  <- unname(object$fixef.means[[k]])
      post_sd    <- apply(draws, 2L, stats::sd)
      n          <- nrow(draws)
      prior_mean <- unname(object$prior$prior_list[[k]]$mu_fixef)
      pval2 <- vapply(seq_len(n_p), function(j) {
        p1 <- mean(draws[, j] < prior_mean[j])
        min(p1, 1 - p1)
      }, numeric(1))
      out <- cbind(
        out,
        fixef.means      = post_mean,
        Post.Sd          = post_sd,
        MC.Error         = post_sd / sqrt(n),
        `Pr(Prior_tail)` = pval2
      )
    }

    rownames(out) <- paste0(k, "::", par)
    out
  })

  do.call(rbind, rows_list)
}

## Per-component tau^2 summary: prior type and plug-in value always; posterior
## mean/SD/quantiles from fixef.dispersion and average candidates per accepted
## Block 2 draw (fixef.iters.mean) for sampled fits.
#' @keywords internal
.lmerb_tau2_summary <- function(object, simulated) {

  ptypes <- object$prior$ptypes
  if (is.null(ptypes)) {
    return(NULL)
  }
  re_names <- object$model_setup$re_coef_names

  tab <- data.frame(
    prior = unname(vapply(re_names, function(k) {
      if (identical(ptypes[[k]], "dIndependent_Normal_Gamma")) "ING" else "dNormal"
    }, character(1L))),
    tau2.plugin = unname(vapply(re_names, function(k) {
      as.numeric(object$prior$prior_list[[k]]$dispersion_fixef)
    }, numeric(1L))),
    row.names = re_names,
    stringsAsFactors = FALSE
  )

  td <- object$fixef.dispersion
  if (simulated && !is.null(td)) {
    tab$Post.Mean <- colMeans(td)[re_names]
    tab$Post.Sd   <- apply(td[, re_names, drop = FALSE], 2L, stats::sd)
    qs <- t(apply(td[, re_names, drop = FALSE], 2L, stats::quantile,
                  probs = c(0.025, 0.5, 0.975)))
    tab$`2.5%`  <- qs[, 1L]
    tab$Median  <- qs[, 2L]
    tab$`97.5%` <- qs[, 3L]
  }
  ## Average envelope candidates per accepted Block 2 draw (1 for dNormal;
  ## ~1/acceptance-rate for ING): sampler-efficiency diagnostic.
  if (simulated && !is.null(object$fixef.iters.mean)) {
    tab$`Cand/draw` <- unname(object$fixef.iters.mean[re_names])
  }

  tab
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

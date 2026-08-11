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
#'   \code{popef} (per-RE-component tables; not printed, available on the
#'   returned object),
#'   \code{ranef_overview}, \code{ranef.iters.mean} (Block~1 envelope candidates
#'   per inner sweep, averaged over groups; printed separately from the overview
#'   table), \code{any_non_normal}, \code{tau2_prior_overview} (per-component
#'   \eqn{\tau^2_k} prior reference: Block~2 \code{pfamily} name (\code{Prior}),
#'   \eqn{1/E[1/\tau^2]}, \eqn{E[\tau^2]}, truncation window,
#'   \code{sqrt(E[tau2])}, and \code{lmer}/\code{glmer} MLE),
#'   \code{sim_method_used} (\code{"DEFAULT"} for the exact-iid engine,
#'   \code{"TWO_BLOCK_GIBBS"} for the two-block Gibbs engine, or \code{NULL}
#'   when \code{simulated} is \code{FALSE}; see \code{\link{lmerb}}'s
#'   \code{sim_method}),
#'   \code{tau2_overview} and \code{tau2_percentiles_overview} (posterior mode,
#'   mean, SD on the variance scale, \code{Mean SD}, and tau^2 quantiles when
#'   simulated), \code{tau2_sd_percentiles_overview} (2.5\%/median/97.5\% of
#'   sqrt(tau^2) draws vs \code{lmer}/\code{glmer SD}), a \code{Residual} row
#'   for observation-level \eqn{\sigma^2} when \code{dispersion_ranef} is
#'   supplied (same columns as the \eqn{\tau^2_k} rows; not used for
#'   per-group \code{gamma_list} priors --- see \code{\link{summary_sigma2}}),
#'   and optionally \code{ranef_groups}.
#' @seealso \code{\link{lmerb}}, \code{\link{glmerb}}, \code{\link{print.lmerb}},
#'   \code{\link{summary_sigma2}}, \code{\link[glmbayes]{summary.glmb}},
#'   \code{\link[glmbayes]{summary.mlmb}}
#' @export
#' @method summary lmerb
summary.lmerb <- function(object, groups = NULL, digits = max(3L, getOption("digits") - 3L), ...) {

  if (!inherits(object, c("lmerb", "glmerb"))) {
    stop("'object' must be an lmerb or glmerb fit.", call. = FALSE)
  }

  re_names  <- object$model_setup$groupef.names
  simulated <- !is.null(object$groupef)
  n_draws   <- if (simulated) nrow(object$popef[[re_names[1L]]]) else NULL
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
    varcor        = VarCorr(mer_fit),
    dispersion    = object$prior$group.dispersion,
    n_obs         = length(object$model_setup$y),
    n_groups      = nlevels(object$model_setup$group),
    group_name    = object$model_setup$group_name,
    fixef_prior_overview = .lmerb_fixef_prior_overview(fixef_parts),
    fixef_overview = .lmerb_fixef_overview(object, simulated = simulated),
    fixef_percentiles_overview = .lmerb_fixef_percentiles_overview(fixef_parts),
    fixef         = fixef_parts,
    ranef_overview = .lmerb_ranef_overview(object, simulated = simulated),
    any_non_normal = isTRUE(object$any_non_normal) ||
      isTRUE(object$prior$any_non_normal),
    tau2_prior_overview       = .lmerb_tau2_prior_overview(object),
    tau2_overview             = .lmerb_tau2_posterior_overview(
      object, simulated = simulated, n_draws = n_draws
    ),
    tau2_percentiles_overview = .lmerb_tau2_percentiles_overview(
      object, simulated = simulated
    ),
    tau2_sd_percentiles_overview = .lmerb_tau2_sd_percentiles_overview(
      object, simulated = simulated
    ),
    ranef.iters.mean = if (simulated) object$groupef.iters.mean else NULL,
    sim_method_used  = object$sim_method_used
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
    cat(sprintf(
      "Bayesian linear mixed model fit  [%d draws, %s]\n",
      x$n, .lmerb_engine_label(x$sim_method_used)))
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

  if (!is.null(x$tau2_prior_overview) && nrow(x$tau2_prior_overview) > 0L) {
    cat("=== Block 2 dispersion (RE variance tau^2_k) ===\n\n")
    cat(sprintf("Prior and %s reference:\n\n", mer_label))
    .lmerb_print_summary_table(x$tau2_prior_overview, digits = digits)
    if (!is.null(x$tau2_overview) && nrow(x$tau2_overview) > 0L) {
      cat("\nOverview:\n")
      stats::printCoefmat(
        round(x$tau2_overview, digits = digits),
        digits = digits,
        quote = FALSE
      )
    }
    if (isTRUE(x$simulated)) {
      if (!is.null(x$tau2_percentiles_overview) &&
          nrow(x$tau2_percentiles_overview) > 0L) {
        cat("\nDistribution percentiles (tau^2):\n\n")
        stats::printCoefmat(
          round(x$tau2_percentiles_overview, digits = digits),
          digits = digits,
          quote = FALSE
        )
      }
      if (!is.null(x$tau2_sd_percentiles_overview) &&
          nrow(x$tau2_sd_percentiles_overview) > 0L) {
        cat(sprintf(
          "\nSD credible interval (sqrt(tau^2) draws; %s SD for reference):\n\n",
          mer_label
        ))
        stats::printCoefmat(
          round(x$tau2_sd_percentiles_overview, digits = digits),
          digits = digits,
          quote = FALSE
        )
      }
    } else {
      cat("\n  (Run with simulate = TRUE for MCMC means, SDs, and percentiles.)\n")
    }
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
  if (isTRUE(x$simulated) && !is.null(x$groupef.iters.mean)) {
    cat(
      "\nMean Block 1 likelihood subgradient candidates per stored draw:",
      formatC(x$groupef.iters.mean, digits = digits, format = "f"),
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
      "Per-group random effects: inspect fit$groupef.mode or fit$groupef,\n",
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
.lmerb_print_summary_table <- function(tab, digits) {
  if (is.null(tab) || nrow(tab) == 0L) {
    return(invisible(tab))
  }
  out <- tab
  num_cols <- vapply(out, is.numeric, logical(1L))
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  print(out, right = TRUE)
  invisible(out)
}

#' Map Block~2 hyper parameter to \code{lmer}/\code{glmer} fixed-effect name
#'
#' Cross-level moderation columns in \code{X_hyper[[slope]]} (e.g. \code{transit_c}
#' moderating \code{log_price_c}) map to the interaction term
#' (\code{transit_c:log_price_c}), not the level-2 main effect \code{transit_c}.
#'
#' @param re_slope_moderation Optional \code{model_setup$re_slope_moderation} data
#'   frame (\code{interaction_col}, \code{moderator}, \code{random_slope}).
#' @keywords internal
.lmerb_lmer_fixef_lookup <- function(
    lmer_fit,
    re_name,
    par_name,
    re_slope_moderation = NULL
) {
  is_glmmtmb <- inherits(lmer_fit, "glmmTMB")
  fe <- if (is_glmmtmb) glmmTMB::fixef(lmer_fit)$cond else fixef(lmer_fit)
  fe_names <- names(fe)

  candidates <- character(0)
  if (par_name == "(Intercept)" && re_name == "(Intercept)") {
    candidates <- "(Intercept)"
  } else if (par_name == "(Intercept)" && re_name != "(Intercept)") {
    candidates <- re_name
  } else if (re_name == "(Intercept)") {
    candidates <- par_name
  } else {
    if (!is.null(re_slope_moderation) && nrow(re_slope_moderation) > 0L) {
      hits <- re_slope_moderation$random_slope == re_name &
        re_slope_moderation$moderator == par_name
      if (any(hits)) {
        candidates <- unique(re_slope_moderation$interaction_col[hits])
      }
    }
    if (!length(candidates)) {
      candidates <- c(
        paste0(par_name, ":", re_name),
        paste0(re_name, ":", par_name),
        par_name
      )
    }
  }

  hit <- candidates[candidates %in% fe_names]
  if (length(hit) == 0L) {
    return(list(estimate = NA_real_, se = NA_real_))
  }

  nm <- hit[1L]
  sm <- tryCatch(summary(lmer_fit), error = function(e) NULL)
  est <- unname(fe[[nm]])
  se  <- NA_real_
  if (!is.null(sm)) {
    coefs <- if (is_glmmtmb) sm$coefficients$cond else sm$coefficients
    if (!is.null(coefs) && nm %in% rownames(coefs)) {
      se <- coefs[nm, "Std. Error"]
    }
  }
  list(estimate = est, se = se)
}

#' @keywords internal
.lmerb_fixef_component_summary <- function(object, k, n_draws, simulated) {

  pl_k   <- object$prior$pop.prior_list[[k]]
  par    <- names(object$popef.mode[[k]])
  q_k    <- length(par)

  prior_mean <- unname(pl_k$mu)
  prior_sd   <- sqrt(diag(pl_k$Sigma))

  mer_ref <- lapply(par, function(nm) {
    .lmerb_lmer_fixef_lookup(
      .lmerb_reference_fit(object),
      k,
      nm,
      re_slope_moderation = object$model_setup$popef.moderation
    )
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

  post_mode <- unname(object$popef.mode[[k]])

  if (!simulated) {
    TAB <- cbind("Post.Mode" = post_mode)
    rownames(TAB) <- par
    return(list(
      coefficients1 = Tab1,
      coefficients  = TAB,
      Percentiles   = NULL
    ))
  }

  draws <- object$popef[[k]]
  post_mean <- unname(object$popef.means[[k]])
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

  re_names <- object$model_setup$groupef.names

  rows_list <- lapply(re_names, function(k) {
    par <- names(object$popef.mode[[k]])
    n_p <- length(par)

    out <- data.frame(
      fixef.mode = unname(object$popef.mode[[k]]),
      stringsAsFactors = FALSE
    )

    if (simulated) {
      draws      <- object$popef[[k]]
      post_mean  <- unname(object$popef.means[[k]])
      post_sd    <- apply(draws, 2L, stats::sd)
      n          <- nrow(draws)
      prior_mean <- unname(object$prior$pop.prior_list[[k]]$mu)
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

## Per-component tau^2 prior reference table.
#' @keywords internal
.lmerb_tau2_prior_overview <- function(object) {

  ptypes <- object$prior$ptypes
  if (is.null(ptypes)) {
    return(NULL)
  }

  re_names  <- object$model_setup$groupef.names
  mer_label <- if (inherits(object, "glmerb")) "glmer" else "lmer"
  mer_vc    <- tryCatch(
    lmebayesCore:::.lmebayes_extract_reference_variance_components(
      .lmerb_reference_fit(object),
      groupef.names = re_names,
      group_name    = object$model_setup$group_name
    ),
    error = function(e) NULL
  )
  vcov_re <- if (!is.null(mer_vc)) mer_vc$Psi else object$model_setup$Psi

  tab <- do.call(rbind, lapply(re_names, function(k) {
    ptype <- ptypes[[k]]
    pf    <- object$prior$pfamily_list[[k]]
    pl    <- if (!is.null(pf)) pf$prior_list else object$prior$pop.prior_list[[k]]
    prior_label <- as.character(ptype)
    mer_tau2 <- if (!is.null(vcov_re) && k %in% names(vcov_re)) {
      unname(vcov_re[[k]])
    } else {
      NA_real_
    }
    mer_sd <- if (is.finite(mer_tau2) && mer_tau2 >= 0) {
      sqrt(mer_tau2)
    } else {
      NA_real_
    }

    if (identical(ptype, "dIndependent_Normal_Gamma")) {
      shape <- as.numeric(pl$shape[1L])
      rate  <- as.numeric(pl$rate[1L])
      inv_E <- if (is.finite(shape) && shape > 0 &&
                    is.finite(rate) && rate > 0) {
        rate / shape
      } else {
        NA_real_
      }
      E_tau2 <- if (is.finite(shape) && shape > 1 &&
                     is.finite(rate) && rate > 0) {
        rate / (shape - 1)
      } else {
        NA_real_
      }
      d_lo <- suppressWarnings(as.numeric(pl$disp_lower))
      d_hi <- suppressWarnings(as.numeric(pl$disp_upper))
      if (!is.finite(d_lo)) d_lo <- NA_real_
      if (!is.finite(d_hi)) d_hi <- NA_real_
    } else {
      d <- as.numeric(object$prior$pop.prior_list[[k]]$dispersion)
      inv_E <- E_tau2 <- d
      d_lo <- d_hi <- NA_real_
    }

    sqrt_E_tau2 <- if (is.finite(E_tau2) && E_tau2 >= 0) {
      sqrt(E_tau2)
    } else {
      NA_real_
    }

    df <- data.frame(
      Prior           = prior_label,
      `1/E[1/tau2]`   = inv_E,
      `E[tau2]`       = E_tau2,
      `sqrt(E[tau2])` = sqrt_E_tau2,
      disp_lower      = d_lo,
      disp_upper      = d_hi,
      check.names     = FALSE,
      stringsAsFactors = FALSE
    )
    df[[mer_label]] <- mer_tau2
    df[[paste0(mer_label, " SD")]] <- mer_sd
    df
  }))

  rownames(tab) <- re_names
  .lmerb_tau2_append_sigma2_row(tab, object, kind = "prior")
}

## Per-group observation-level sigma^2 when dispersion_ranef is a gamma_list.
#' @keywords internal
.lmerb_sigma2_gamma_list_prior_overview <- function(object) {
  if (!identical(object$prior$dispersion_mode, "gamma_list")) {
    return(NULL)
  }
  pl <- object$prior$dispersion_prior_list
  if (is.null(pl) || is.null(pl$shape_group)) {
    return(NULL)
  }
  grp <- names(pl$shape_group)
  if (!length(grp)) {
    return(NULL)
  }
  mer <- .lmerb_sigma2_mer_reference(object)

  ## Number of random-effect coefficients (p) shared by every group's dGamma()
  ## envelope; used to back out the implied per-group effective prior sample
  ## size n_prior under the Prior_Setup() default k = 1 calibration
  ## (shape_ING = (n_prior + k + p)/2), i.e. the same convention used by
  ## lmebayesCore:::.ing_n_prior_from_shape().
  p_re <- length(object$model_setup$groupef.names)
  n_data_group <- table(object$model_setup$group)

  tab <- do.call(rbind, lapply(grp, function(g) {
    sh <- as.numeric(pl$shape_group[[g]])
    rt <- as.numeric(pl$rate_group[[g]])
    n_prior_g <- if (is.finite(sh)) 2 * sh - 1 - p_re else NA_real_
    n_data_g  <- if (g %in% names(n_data_group)) {
      as.numeric(n_data_group[[g]])
    } else {
      NA_real_
    }
    inv_E <- if (is.finite(sh) && sh > 0 && is.finite(rt) && rt > 0) {
      rt / sh
    } else {
      NA_real_
    }
    E_sigma2 <- if (is.finite(sh) && sh > 1 && is.finite(rt) && rt > 0) {
      rt / (sh - 1)
    } else {
      NA_real_
    }
    sqrt_E <- if (is.finite(E_sigma2) && E_sigma2 >= 0) {
      sqrt(E_sigma2)
    } else {
      NA_real_
    }
    d_lo <- suppressWarnings(as.numeric(pl$disp_lower_group[[g]]))
    d_hi <- suppressWarnings(as.numeric(pl$disp_upper_group[[g]]))
    if (!is.finite(d_lo)) d_lo <- NA_real_
    if (!is.finite(d_hi)) d_hi <- NA_real_

    wdiag <- object$prior$window_diagnostics
    wrow <- if (!is.null(wdiag) && g %in% wdiag$group) {
      wdiag[wdiag$group == g, , drop = FALSE]
    } else {
      NULL
    }

    df <- data.frame(
      Prior             = "dGamma",
      n_prior           = n_prior_g,
      `1/E[1/sigma2]`   = inv_E,
      `E[sigma2]`       = E_sigma2,
      `sqrt(E[sigma2])` = sqrt_E,
      disp_lower        = d_lo,
      disp_upper        = d_hi,
      n_data            = n_data_g,
      check.names       = FALSE,
      stringsAsFactors  = FALSE
    )
    if (!is.null(wrow) && nrow(wrow) == 1L) {
      df$blup_infl <- wrow$blup_infl
      df$R_lo <- wrow$R_lo
      df$R_hi <- wrow$R_hi
      df$asymmetric_window <- wrow$asymmetric_window
    }
    df[[mer$mer_label]] <- mer$mer_sigma2
    df[[paste0(mer$mer_label, " SD")]] <- mer$mer_sd
    df
  }))
  rownames(tab) <- grp
  tab
}

#' @keywords internal
.lmerb_sigma2_gamma_list_posterior_overview <- function(
    object,
    simulated = FALSE,
    n_draws = NULL
) {
  prior_tab <- .lmerb_sigma2_gamma_list_prior_overview(object)
  if (is.null(prior_tab) || nrow(prior_tab) == 0L) {
    return(NULL)
  }
  grp <- rownames(prior_tab)
  E_prior <- prior_tab[["E[sigma2]"]]
  n_total <- prior_tab[["n_prior"]] + prior_tab[["n_data"]]

  post_mode <- if (!is.null(object$group.dispersion.mode)) {
    sm <- object$group.dispersion.mode
    if (is.matrix(sm)) {
      vapply(grp, function(g) {
        j <- match(g, colnames(sm))
        if (is.na(j)) NA_real_ else as.numeric(sm[1L, j])
      }, numeric(1L))
    } else {
      v <- as.numeric(sm)
      if (length(v) == length(grp)) {
        v[match(grp, names(v))]
      } else if (length(v) == 1L) {
        rep(v, length(grp))
      } else {
        rep(NA_real_, length(grp))
      }
    }
  } else {
    rep(NA_real_, length(grp))
  }

  if (!isTRUE(simulated)) {
    out <- cbind(`n_total` = n_total, `Post.Mode` = post_mode)
    rownames(out) <- grp
    return(out)
  }

  sigma2 <- object$group.dispersion
  if (is.null(sigma2) || !is.matrix(sigma2)) {
    out <- cbind(`n_total` = n_total, `Post.Mode` = post_mode)
    rownames(out) <- grp
    return(out)
  }

  tab <- do.call(rbind, lapply(seq_along(grp), function(i) {
    g <- grp[i]
    j <- match(g, colnames(sigma2))
    draws <- if (!is.na(j)) sigma2[, j] else NA_real_
    post_mean <- mean(draws)
    post_sd   <- stats::sd(draws)
    mean_sd   <- mean(sqrt(draws))
    mc_err    <- if (!is.null(n_draws) && n_draws > 0L) {
      post_sd / sqrt(n_draws)
    } else {
      NA_real_
    }
    pval2 <- if (is.finite(E_prior[i])) {
      p1 <- mean(draws < E_prior[i])
      min(p1, 1 - p1)
    } else {
      NA_real_
    }
    c(
      `n_total`        = n_total[i],
      `Post.Mode`      = post_mode[i],
      `Post.Mean`      = post_mean,
      `Post.Sd`        = post_sd,
      `MC Error`       = mc_err,
      `Mean SD`        = mean_sd,
      `Pr(Prior_tail)` = pval2
    )
  }))
  rownames(tab) <- grp

  cand <- if (!is.null(object$group.dispersion.iters.mean)) {
    v <- object$group.dispersion.iters.mean[grp]
    if (is.null(names(v))) {
      stats::setNames(v, grp)
    } else {
      v
    }
  } else if (!is.null(object$groupef.iters.mean)) {
    rep(unname(object$groupef.iters.mean), length(grp))
  } else {
    NULL
  }
  if (!is.null(cand)) {
    tab <- cbind(tab, `Cand/draw` = unname(cand))
  }

  tab
}

## Per-group sigma^2 posterior percentiles, for comparison against
## disp_lower/disp_upper (see .lmerb_sigma2_gamma_list_prior_overview()) --
## mirrors .lmerb_tau2_percentiles_overview()'s Block~2 layout, so it is
## easy to see whether draws pile up against the truncation bounds.
#' @keywords internal
.lmerb_sigma2_gamma_list_percentiles_overview <- function(object, simulated) {
  if (!isTRUE(simulated)) {
    return(NULL)
  }
  sigma2 <- object$group.dispersion
  if (is.null(sigma2) || !is.matrix(sigma2)) {
    return(NULL)
  }
  grp <- colnames(sigma2)
  if (is.null(grp) || !length(grp)) {
    return(NULL)
  }

  percentiles <- t(apply(sigma2, 2L, stats::quantile,
    probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
  ))
  tab <- cbind(
    `1.0%`   = percentiles[, 1L],
    `2.5%`   = percentiles[, 2L],
    `5.0%`   = percentiles[, 3L],
    Median   = percentiles[, 4L],
    `95.0%`  = percentiles[, 5L],
    `97.5%`  = percentiles[, 6L],
    `99.0%`  = percentiles[, 7L]
  )
  rownames(tab) <- grp
  tab
}

## Observation-level sigma^2: append a Residual row to tau^2 summary tables.
#' @keywords internal
.lmerb_sigma2_summary_enabled <- function(object) {
  mode <- object$prior$dispersion_mode
  # "gamma_list" (per-group dGamma() priors) draws sigma2 as an n x J
  # matrix, and "fixed_vector" (per-group known constants) is a length-J
  # vector; the pooled prior-vs-posterior Residual row below only makes
  # sense for a single pooled value, so both per-group modes are routed to
  # summary_sigma2() instead (see .lmerb_sigma2_fixed_vector_overview() and
  # .lmerb_sigma2_gamma_list_prior_overview()).
  !is.null(mode) && !identical(mode, "none") &&
    !mode %in% c("gamma_list", "fixed_vector")
}

## Per-group observation-level sigma^2 when dispersion_ranef is a fixed,
## known numeric vector (mode = "fixed_vector"): no prior to calibrate and
## nothing sampled, so this is just the constant values themselves next to
## the lmer/glmer residual variance for reference.
#' @keywords internal
.lmerb_sigma2_fixed_vector_overview <- function(object) {
  if (!identical(object$prior$dispersion_mode, "fixed_vector")) {
    return(NULL)
  }
  sigma2 <- object$prior$group.dispersion
  if (is.null(sigma2) || !length(sigma2)) {
    return(NULL)
  }
  grp <- names(sigma2)
  if (is.null(grp) || !length(grp)) {
    grp <- levels(object$model_setup$group)
  }
  mer <- .lmerb_sigma2_mer_reference(object)

  df <- data.frame(
    `Fixed sigma2` = as.numeric(sigma2),
    check.names    = FALSE,
    stringsAsFactors = FALSE
  )
  df[[mer$mer_label]] <- mer$mer_sigma2
  df[[paste0(mer$mer_label, " SD")]] <- mer$mer_sd
  rownames(df) <- grp
  df
}

#' Pooled \code{lmer}/\code{glmer} residual variance, for comparison against
#' per-group dispersion tables (\code{gamma_list} / \code{fixed_vector}).
#'
#' Deliberately uses the plain pooled reference fit
#' (\code{object$lmer}/\code{object$glmer}), \strong{not}
#' \code{\link{.lmerb_reference_fit}}: when \code{dispformula} requests
#' per-group dispersion, the whole point of these tables is comparing the
#' per-group values against a single pooled baseline, and a \code{glmmTMB}
#' fit with a per-group \code{dispformula} has no single residual variance
#' to report (\code{extract_glmmtmb_variance_components()} always returns
#' \code{residual_var = NA} for it).
#' @keywords internal
#' @noRd
.lmerb_sigma2_mer_reference <- function(object) {
  re_names   <- object$model_setup$groupef.names
  mer_label  <- if (inherits(object, "glmerb")) "glmer" else "lmer"
  pooled_fit <- if (inherits(object, "glmerb")) object$glmer else object$lmer
  mer_vc    <- tryCatch(
    lmebayesCore:::extract_mer_variance_components(
      pooled_fit,
      groupef.names = re_names
    ),
    error = function(e) NULL
  )
  resid_var <- if (!is.null(mer_vc)) {
    mer_vc$dispersion
  } else {
    object$model_setup$dispersion
  }
  mer_sd <- if (is.finite(resid_var) && resid_var >= 0) {
    sqrt(resid_var)
  } else {
    NA_real_
  }
  list(
    mer_label  = mer_label,
    mer_sigma2 = resid_var,
    mer_sd     = mer_sd
  )
}

#' @keywords internal
.lmerb_sigma2_prior_params <- function(object) {
  mode <- object$prior$dispersion_mode
  pl   <- object$prior$dispersion_prior_list
  mer  <- .lmerb_sigma2_mer_reference(object)

  if (identical(mode, "gamma") && !is.null(pl)) {
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    inv_E <- if (is.finite(shape) && shape > 0 &&
                  is.finite(rate) && rate > 0) {
      rate / shape
    } else {
      NA_real_
    }
    E_sigma2 <- if (is.finite(shape) && shape > 1 &&
                       is.finite(rate) && rate > 0) {
      rate / (shape - 1)
    } else {
      NA_real_
    }
    d_lo <- suppressWarnings(as.numeric(pl$disp_lower))
    d_hi <- suppressWarnings(as.numeric(pl$disp_upper))
    if (!is.finite(d_lo)) d_lo <- NA_real_
    if (!is.finite(d_hi)) d_hi <- NA_real_
    prior_label <- "dGamma"
  } else if (identical(mode, "fixed")) {
    d <- as.numeric(object$prior$group.dispersion)
    inv_E <- E_sigma2 <- d
    d_lo <- d_hi <- NA_real_
    prior_label <- "fixed"
  } else {
    return(NULL)
  }

  sqrt_E <- if (is.finite(E_sigma2) && E_sigma2 >= 0) {
    sqrt(E_sigma2)
  } else {
    NA_real_
  }

  list(
    prior_label = prior_label,
    inv_E       = inv_E,
    E_sigma2    = E_sigma2,
    sqrt_E      = sqrt_E,
    disp_lower  = d_lo,
    disp_upper  = d_hi,
    mer_sigma2  = mer$mer_sigma2,
    mer_sd      = mer$mer_sd,
    mer_label   = mer$mer_label
  )
}

#' @keywords internal
.lmerb_tau2_append_sigma2_row <- function(
    tab,
    object,
    kind = c("prior", "overview", "percentiles", "sd_percentiles"),
    simulated = FALSE,
    n_draws = NULL
) {
  kind <- match.arg(kind)
  if (!.lmerb_sigma2_summary_enabled(object)) {
    return(tab)
  }
  if (is.null(tab) || nrow(tab) == 0L) {
    return(tab)
  }

  params <- .lmerb_sigma2_prior_params(object)
  if (is.null(params)) {
    return(tab)
  }

  if (identical(kind, "prior")) {
    df <- data.frame(
      Prior           = params$prior_label,
      `1/E[1/tau2]`   = params$inv_E,
      `E[tau2]`       = params$E_sigma2,
      `sqrt(E[tau2])` = params$sqrt_E,
      disp_lower      = params$disp_lower,
      disp_upper      = params$disp_upper,
      check.names     = FALSE,
      stringsAsFactors = FALSE
    )
    df[[params$mer_label]] <- params$mer_sigma2
    df[[paste0(params$mer_label, " SD")]] <- params$mer_sd
    rownames(df) <- "Residual"
    return(rbind(tab, df))
  }

  if (identical(kind, "overview")) {
    post_mode <- if (identical(object$prior$dispersion_mode, "fixed")) {
      as.numeric(object$prior$group.dispersion)
    } else {
      params$E_sigma2
    }
    if (!simulated) {
      out <- cbind(`Post.Mode` = post_mode)
      rownames(out) <- "Residual"
      return(rbind(tab, out))
    }

    sigma2 <- object$group.dispersion
    if (is.null(sigma2)) {
      out <- cbind(`Post.Mode` = post_mode)
      rownames(out) <- "Residual"
      return(rbind(tab, out))
    }
    sigma2 <- as.numeric(sigma2)
    post_mean <- if (!is.null(object$group.dispersion.mean)) {
      as.numeric(object$group.dispersion.mean)
    } else {
      mean(sigma2)
    }
    post_sd <- stats::sd(sigma2)
    mean_sd <- mean(sqrt(sigma2))
    mc_err <- if (!is.null(n_draws) && n_draws > 0L) {
      post_sd / sqrt(n_draws)
    } else {
      NA_real_
    }
    pval2 <- if (identical(object$prior$dispersion_mode, "gamma") &&
                  is.finite(params$E_sigma2)) {
      p1 <- mean(sigma2 < params$E_sigma2)
      min(p1, 1 - p1)
    } else {
      NA_real_
    }
    out <- cbind(
      `Post.Mode`      = post_mode,
      `Post.Mean`      = post_mean,
      `Post.Sd`        = post_sd,
      `MC Error`       = mc_err,
      `Mean SD`        = mean_sd,
      `Pr(Prior_tail)` = pval2
    )
    if ("Cand/draw" %in% colnames(tab)) {
      ## 'tab' already has this column (added upstream whenever
      ## fixef.iters.mean is non-NULL); rbind() requires the sigma2 row to
      ## have it too, even when the Block~1 rate (ranef.iters.mean) is NA
      ## (e.g. dispersion_ranef is fixed/known -- no Block~1 envelope).
      cand_val <- if (!is.null(object$groupef.iters.mean)) {
        unname(object$groupef.iters.mean)
      } else {
        NA_real_
      }
      out <- cbind(out, `Cand/draw` = cand_val)
    }
    rownames(out) <- "Residual"
    return(rbind(tab, out))
  }

  if (identical(kind, "percentiles")) {
    if (!simulated) {
      return(tab)
    }
    sigma2 <- object$group.dispersion
    if (is.null(sigma2)) {
      return(tab)
    }
    sigma2 <- as.numeric(sigma2)
    pct <- stats::quantile(
      sigma2,
      probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
    )
    out <- cbind(
      `1.0%`   = pct[[1L]],
      `2.5%`   = pct[[2L]],
      `5.0%`   = pct[[3L]],
      Median   = pct[[4L]],
      `95.0%`  = pct[[5L]],
      `97.5%`  = pct[[6L]],
      `99.0%`  = pct[[7L]]
    )
    rownames(out) <- "Residual"
    return(rbind(tab, out))
  }

  if (identical(kind, "sd_percentiles")) {
    if (!simulated) {
      return(tab)
    }
    sigma2 <- object$group.dispersion
    if (is.null(sigma2)) {
      return(tab)
    }
    sigma2 <- as.numeric(sigma2)
    sd_draws <- sqrt(sigma2)
    pct <- stats::quantile(sd_draws, probs = c(0.025, 0.5, 0.975))
    mer_sd_col <- paste0(params$mer_label, " SD")
    out <- cbind(
      `2.5%`  = pct[[1L]],
      Median  = pct[[2L]],
      `97.5%` = pct[[3L]],
      params$mer_sd
    )
    colnames(out)[ncol(out)] <- mer_sd_col
    rownames(out) <- "Residual"
    return(rbind(tab, out))
  }

  tab
}

## Per-component tau^2 posterior overview (mode at plug-in / fixed value;
## MCMC mean, SD, tail probability vs E[tau2], envelope candidates).
#' @keywords internal
.lmerb_tau2_posterior_overview <- function(object, simulated, n_draws) {

  ptypes <- object$prior$ptypes
  if (is.null(ptypes)) {
    return(NULL)
  }

  re_names <- object$model_setup$groupef.names
  prior_tab <- .lmerb_tau2_prior_overview(object)
  if (is.null(prior_tab)) {
    return(NULL)
  }

  post_mode <- vapply(re_names, function(k) {
    as.numeric(object$prior$pop.prior_list[[k]]$dispersion)
  }, numeric(1))

  if (!simulated) {
    out <- cbind(`Post.Mode` = post_mode)
    rownames(out) <- re_names
    return(.lmerb_tau2_append_sigma2_row(
      out, object, kind = "overview", simulated = simulated, n_draws = n_draws
    ))
  }

  td <- object$popef.dispersion
  if (is.null(td)) {
    out <- cbind(`Post.Mode` = post_mode)
    rownames(out) <- re_names
    return(.lmerb_tau2_append_sigma2_row(
      out, object, kind = "overview", simulated = simulated, n_draws = n_draws
    ))
  }

  post_mean <- colMeans(td)[re_names]
  post_sd   <- apply(td[, re_names, drop = FALSE], 2L, stats::sd)
  mean_sd   <- vapply(re_names, function(k) {
    mean(sqrt(td[, k]))
  }, numeric(1))
  mc_err    <- if (!is.null(n_draws) && n_draws > 0L) {
    post_sd / sqrt(n_draws)
  } else {
    rep(NA_real_, length(re_names))
  }

  E_prior <- prior_tab[re_names, "E[tau2]", drop = TRUE]
  pval2 <- vapply(seq_along(re_names), function(j) {
    k <- re_names[j]
    if (identical(ptypes[[k]], "dNormal")) {
      return(NA_real_)
    }
    p1 <- mean(td[, k] < E_prior[j])
    min(p1, 1 - p1)
  }, numeric(1))

  out <- cbind(
    `Post.Mode`      = post_mode,
    `Post.Mean`      = post_mean,
    `Post.Sd`        = post_sd,
    `MC Error`       = mc_err,
    `Mean SD`        = mean_sd,
    `Pr(Prior_tail)` = pval2
  )

  if (!is.null(object$popef.iters.mean)) {
    out <- cbind(out, `Cand/draw` = unname(object$popef.iters.mean[re_names]))
  }

  rownames(out) <- re_names
  .lmerb_tau2_append_sigma2_row(
    out, object, kind = "overview", simulated = simulated, n_draws = n_draws
  )
}

## Per-component tau^2 posterior percentiles from fixef.dispersion draws.
#' @keywords internal
.lmerb_tau2_percentiles_overview <- function(object, simulated) {

  if (!simulated) {
    return(NULL)
  }

  ptypes <- object$prior$ptypes
  if (is.null(ptypes)) {
    return(NULL)
  }

  re_names <- object$model_setup$groupef.names
  td <- object$popef.dispersion
  if (is.null(td)) {
    return(NULL)
  }

  percentiles <- t(apply(td[, re_names, drop = FALSE], 2L, stats::quantile,
    probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
  ))
  tab <- cbind(
    `1.0%`   = percentiles[, 1L],
    `2.5%`   = percentiles[, 2L],
    `5.0%`   = percentiles[, 3L],
    Median   = percentiles[, 4L],
    `95.0%`  = percentiles[, 5L],
    `97.5%`  = percentiles[, 6L],
    `99.0%`  = percentiles[, 7L]
  )
  rownames(tab) <- re_names
  .lmerb_tau2_append_sigma2_row(
    tab, object, kind = "percentiles", simulated = simulated
  )
}

## Per-component SD (sqrt(tau^2)) posterior percentiles from fixef.dispersion draws.
#' @keywords internal
.lmerb_tau2_sd_percentiles_overview <- function(object, simulated) {

  if (!simulated) {
    return(NULL)
  }

  ptypes <- object$prior$ptypes
  if (is.null(ptypes)) {
    return(NULL)
  }

  re_names <- object$model_setup$groupef.names
  td <- object$popef.dispersion
  if (is.null(td)) {
    return(NULL)
  }

  prior_tab <- .lmerb_tau2_prior_overview(object)
  mer_label <- if (inherits(object, "glmerb")) "glmer" else "lmer"
  mer_sd_col <- paste0(mer_label, " SD")
  mer_sd <- if (!is.null(prior_tab) && mer_sd_col %in% colnames(prior_tab)) {
    prior_tab[re_names, mer_sd_col, drop = TRUE]
  } else {
    rep(NA_real_, length(re_names))
  }

  sd_draws <- sqrt(td[, re_names, drop = FALSE])
  percentiles <- t(apply(sd_draws, 2L, stats::quantile,
    probs = c(0.025, 0.5, 0.975)
  ))
  tab <- cbind(
    `2.5%`  = percentiles[, 1L],
    Median  = percentiles[, 2L],
    `97.5%` = percentiles[, 3L]
  )
  tab <- cbind(tab, mer_sd)
  colnames(tab)[ncol(tab)] <- mer_sd_col
  rownames(tab) <- re_names
  .lmerb_tau2_append_sigma2_row(
    tab, object, kind = "sd_percentiles", simulated = simulated
  )
}

#' @keywords internal
.lmerb_ranef_overview <- function(object, simulated) {

  re_names <- object$model_setup$groupef.names
  b_mode   <- object$groupef.mode

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
        object$groupef[[k]],
        object$groupef[[grp_col]],
        mean
      ))
    }, numeric(1))
    overview <- cbind(overview, MCMC.mean = mcmc_means)
  }

  overview
}

#' @keywords internal
.lmerb_ranef_groups_detail <- function(object, groups, simulated) {

  re_names <- object$model_setup$groupef.names
  grp_col  <- object$model_setup$group_name
  groups   <- as.character(groups)

  rows <- lapply(groups, function(lev) {
    mode_vals <- object$groupef.mode[lev, re_names, drop = TRUE]
    out <- data.frame(
      group = lev,
      t(mode_vals),
      check.names = FALSE
    )
    names(out) <- c("group", re_names)

    if (simulated) {
      idx <- object$groupef[[grp_col]] == lev
      for (k in re_names) {
        out[[paste0(k, ".mean")]] <- mean(object$groupef[idx, k])
        out[[paste0(k, ".sd")]]   <- stats::sd(object$groupef[idx, k])
      }
    }
    out
  })

  do.call(rbind, rows)
}

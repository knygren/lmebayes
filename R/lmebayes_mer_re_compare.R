#' Anchor for lmebayes-aligned mer_full random effects
#'
#' @param re_names Character vector of RE component names.
#' @param mer_fixef Named numeric vector from \code{lme4::fixef(mer_fit)}.
#' @param fixef_mode Optional list of ICM Block~2 modes (fallback for pure RE terms).
#' @return Named numeric vector, one anchor per RE component.
#' @keywords internal
#' @noRd
.mer_coef_anchor <- function(re_names, mer_fixef, fixef_mode = NULL) {
  vapply(re_names, function(k) {
    if (k == "(Intercept)") {
      unname(mer_fixef["(Intercept)"])
    } else if (k %in% names(mer_fixef)) {
      unname(mer_fixef[k])
    } else if (!is.null(fixef_mode)) {
      unname(fixef_mode[[k]]["(Intercept)"])
    } else {
      0
    }
  }, numeric(1L))
}

#' lmebayes-aligned MER random effects at each grouping level
#'
#' \code{mer_full = mu_all + (coef(mer) - anchor)} using the per-group
#' population expectation from \code{.lmerb_popef_mu()}.
#'
#' @param fit An \code{lmerb} or \code{glmerb} object with \code{$lmer} or \code{$glmer}.
#' @return Numeric matrix \code{J x p_re} with rownames = group levels.
#' @keywords internal
#' @noRd
.mer_re_reference_full <- function(fit) {
  if (!inherits(fit, c("lmerb", "glmerb"))) {
    stop("fit must be an lmerb or glmerb object.", call. = FALSE)
  }
  ## .lmerb_reference_fit() prefers fit$glmmTMB (per-group dispersion
  ## reference) over fit$lmer/fit$glmer when present, so comparisons match
  ## whichever fit actually calibrated the priors; .lmebayes_reference_coef()
  ## / .lmebayes_reference_fixef() dispatch coef()/fixef() uniformly for
  ## either fit type (glmmTMB's are nested under $cond).
  mer_fit   <- .lmerb_reference_fit(fit)
  re_names <- fit$model_setup$groupef.names
  grp_col  <- fit$model_setup$group_name
  mu_mat   <- .lmerb_popef_mu(fit)
  if (is.null(mu_mat)) {
    stop("Population expectation is unavailable for this fit.", call. = FALSE)
  }

  coef_raw <- as.data.frame(
    lmebayesCore:::.lmebayes_reference_coef(mer_fit)[[grp_col]][, re_names, drop = FALSE]
  )
  grp_levs <- rownames(coef_raw)
  fe_mer   <- lmebayesCore:::.lmebayes_reference_fixef(mer_fit)
  anchor   <- .mer_coef_anchor(re_names, fe_mer, .lmerb_popef_mode(fit))

  mer_full <- matrix(
    NA_real_,
    nrow = length(grp_levs),
    ncol = length(re_names),
    dimnames = list(grp_levs, re_names)
  )
  for (j in seq_len(nrow(coef_raw))) {
    lev <- grp_levs[j]
    for (k in re_names) {
      mu_k <- mu_mat[k, lev, drop = TRUE]
      mer_full[lev, k] <- mu_k + (unname(coef_raw[[k]][j]) - anchor[k])
    }
  }
  mer_full
}

#' Chain-mean random effects by grouping level
#'
#' @param fit An \code{lmerb} or \code{glmerb} object with stored \code{$groupef}.
#' @return Numeric matrix \code{J x p_re}, or \code{NULL} when no draws.
#' @keywords internal
#' @noRd
.bayes_ranef_chain_means <- function(fit) {
  draws <- .lmerb_groupef_draws(fit)
  if (is.null(draws)) {
    return(NULL)
  }
  re_names <- fit$model_setup$groupef.names
  grp_col  <- fit$model_setup$group_name
  re_draws_mean <- tapply(
    seq_len(nrow(draws)),
    draws[[grp_col]],
    function(idx) colMeans(draws[idx, re_names, drop = FALSE]),
    simplify = FALSE
  )
  grp_levs <- rownames(.lmerb_groupef_mode(fit))
  if (is.null(grp_levs)) {
    grp_levs <- levels(factor(draws[[grp_col]]))
  }
  if (!setequal(names(re_draws_mean), grp_levs)) {
    stop(
      "Group levels in groupef draws do not match the fitted group levels.",
      call. = FALSE
    )
  }
  bayes <- do.call(rbind, re_draws_mean[grp_levs])
  rownames(bayes) <- grp_levs
  colnames(bayes) <- re_names
  bayes
}

#' @keywords internal
#' @noRd
.mer_re_rename_cols <- function(df, re_names, suffix) {
  idx <- match(re_names, names(df))
  if (anyNA(idx)) {
    stop(
      "missing RE columns: ",
      paste(re_names[is.na(idx)], collapse = ", "),
      call. = FALSE
    )
  }
  names(df)[idx] <- paste0(re_names, suffix)
  df
}

#' Print side-by-side MER reference vs Bayesian chain-mean random effects
#'
#' Compares \code{mer_full = mu_all + (coef(mer) - anchor)} to the mean of
#' \code{fit$groupef} within each grouping level. Requires
#' \code{simulate = TRUE} (stored draws).
#'
#' @param fit An \code{lmerb} or \code{glmerb} object.
#' @param detail \code{"compact"} (default) or \code{"full"} (Ex_12-style tables).
#' @param digits Decimal places for printed numerics.
#' @return Invisibly, a list with \code{mer_full}, \code{bayes}, and \code{anchor}.
#' @keywords internal
#' @noRd
print_mer_bayes_re_compare <- function(
    fit,
    detail = c("compact", "full"),
    digits = 4L
) {
  detail <- match.arg(detail)
  if (!inherits(fit, c("lmerb", "glmerb"))) {
    stop("fit must be an lmerb or glmerb object.", call. = FALSE)
  }

  is_lmerb   <- inherits(fit, "lmerb")
  mer_label  <- if (is_lmerb) "lmer" else "glmer"
  bayes_label <- if (is_lmerb) "lmerb" else "glmerb"
  re_names   <- fit$model_setup$groupef.names
  grp_col    <- fit$model_setup$group_name
  popef_draws <- .lmerb_popef_draws(fit)
  n_draws    <- if (!is.null(popef_draws) && length(re_names)) {
    nrow(popef_draws[[re_names[1L]]])
  } else {
    NA_integer_
  }

  bayes <- .bayes_ranef_chain_means(fit)
  if (is.null(bayes)) {
    cat(
      "Random-effects comparison requires simulate = TRUE ",
      "(stored fit$groupef).\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  mer_full <- .mer_re_reference_full(fit)
  mer_fit  <- .lmerb_reference_fit(fit)
  fe_mer   <- lmebayesCore:::.lmebayes_reference_fixef(mer_fit)
  anchor   <- .mer_coef_anchor(re_names, fe_mer, .lmerb_popef_mode(fit))
  grp_levs <- rownames(mer_full)

  cat(sprintf(
    "  %s: colMeans(fit$groupef) within each %s (n = %s chains)\n",
    bayes_label, grp_col,
    if (is.na(n_draws)) "?" else format(n_draws, big.mark = ",")
  ))
  cat(sprintf(
    "  %s:  mer_full = mu_all + (coef(%s) - anchor)\n\n",
    mer_label, mer_label
  ))

  if (detail == "compact") {
    cat(sprintf("  Mean across %s levels:\n", grp_col))
    avg_re <- data.frame(
      term = re_names,
      mer_full = vapply(re_names, function(k) mean(mer_full[, k]), numeric(1)),
      bayes = vapply(re_names, function(k) mean(bayes[, k]), numeric(1)),
      row.names = NULL
    )
    names(avg_re)[2L:3L] <- c(
      paste0(mer_label, "_full"),
      bayes_label
    )
    avg_re[-1L] <- lapply(avg_re[-1L], function(x) round(x, digits))
    print(avg_re)

    cat(sprintf(
      "\n  Per-%s comparison (%s_full vs %s chain mean):\n",
      grp_col, mer_label, bayes_label
    ))
    for (k in re_names) {
      cat(sprintf("\n  --- %s ---\n", k))
      cmp_k <- cbind(
        mer_full = mer_full[, k, drop = TRUE],
        bayes = bayes[, k, drop = TRUE]
      )
      colnames(cmp_k) <- c(
        paste0(mer_label, "_full"),
        bayes_label
      )
      print(round(cmp_k, digits))
    }
  } else {
    coef_raw_df <- as.data.frame(
      lmebayesCore:::.lmebayes_reference_coef(mer_fit)[[grp_col]][, re_names, drop = FALSE]
    )
    mu_mat      <- .lmerb_popef_mu(fit)

    mer_by_level <- coef_raw_df
    mer_by_level[[grp_col]] <- factor(rownames(mer_by_level), levels = grp_levs)
    rownames(mer_by_level) <- NULL
    mer_by_level <- .mer_re_rename_cols(mer_by_level, re_names, paste0("_", mer_label))

    mu_by_level <- as.data.frame(t(mu_mat), stringsAsFactors = FALSE)
    mu_by_level[[grp_col]] <- factor(rownames(mu_by_level), levels = grp_levs)
    rownames(mu_by_level) <- NULL
    mu_by_level <- .mer_re_rename_cols(mu_by_level, re_names, "_mu_all")

    mer_full_df <- as.data.frame(mer_full, stringsAsFactors = FALSE)
    mer_full_by_level <- mer_full_df
    mer_full_by_level[[grp_col]] <- factor(grp_levs, levels = grp_levs)
    rownames(mer_full_by_level) <- NULL
    mer_full_by_level <- .mer_re_rename_cols(
      mer_full_by_level, re_names, paste0("_", mer_label, "_full")
    )

    bayes_df <- as.data.frame(bayes, stringsAsFactors = FALSE)
    bayes_by_level <- bayes_df
    bayes_by_level[[grp_col]] <- factor(grp_levs, levels = grp_levs)
    rownames(bayes_by_level) <- NULL
    bayes_by_level <- .mer_re_rename_cols(
      bayes_by_level, re_names, paste0("_", bayes_label)
    )

    level_means <- merge(mer_by_level, mer_full_by_level, by = grp_col, sort = FALSE)
    level_means <- merge(level_means, mu_by_level, by = grp_col, sort = FALSE)
    level_means <- merge(level_means, bayes_by_level, by = grp_col, sort = FALSE)
    level_means <- level_means[order(level_means[[grp_col]]), , drop = FALSE]

    mer_raw_col   <- paste0(re_names, "_", mer_label)
    mer_full_col  <- paste0(re_names, "_", mer_label, "_full")
    mu_all_col    <- paste0(re_names, "_mu_all")
    bayes_col     <- paste0(re_names, "_", bayes_label)

    cat(sprintf(
      "\nMean coefficient across %s levels (raw coef vs mu_all vs mer_full vs %s):\n",
      grp_col, bayes_label
    ))
    avg_row <- data.frame(
      term = re_names,
      mer_raw = vapply(re_names, function(nm) {
        mean(level_means[[paste0(nm, "_", mer_label)]])
      }, numeric(1L)),
      mu_all = vapply(re_names, function(nm) {
        mean(level_means[[paste0(nm, "_mu_all")]])
      }, numeric(1L)),
      mer_full = vapply(re_names, function(nm) {
        mean(level_means[[paste0(nm, "_", mer_label, "_full")]])
      }, numeric(1L)),
      bayes = vapply(re_names, function(nm) {
        mean(level_means[[paste0(nm, "_", bayes_label)]])
      }, numeric(1L)),
      row.names = NULL
    )
    names(avg_row)[-1L] <- c(
      paste0(mer_label, "_raw"),
      "mu_all",
      paste0(mer_label, "_full"),
      bayes_label
    )
    avg_row[-1L] <- lapply(avg_row[-1L], function(x) round(x, digits))
    print(avg_row)

    show_cols <- c(grp_col, mer_raw_col, mer_full_col, mu_all_col, bayes_col)
    cat(sprintf(
      "\nFactor-level means: %s (raw coef), mer_full, mu_all, %s:\n",
      mer_label, bayes_label
    ))
    out_means <- level_means[, show_cols, drop = FALSE]
    num_cols  <- setdiff(show_cols, grp_col)
    out_means[num_cols] <- lapply(out_means[num_cols], function(x) {
      round(as.numeric(x), digits)
    })
    print(out_means)

    level_long <- do.call(rbind, lapply(re_names, function(nm) {
      mer_full_v <- level_means[[paste0(nm, "_", mer_label, "_full")]]
      mu_all_v   <- level_means[[paste0(nm, "_mu_all")]]
      bayes_v    <- level_means[[paste0(nm, "_", bayes_label)]]
      data.frame(
        level     = level_means[[grp_col]],
        term      = nm,
        mer_raw   = level_means[[paste0(nm, "_", mer_label)]],
        mer_full  = mer_full_v,
        mu_all    = mu_all_v,
        bayes     = bayes_v,
        u_mer     = mer_full_v - mu_all_v,
        u_bayes   = bayes_v - mu_all_v,
        diff_bf   = bayes_v - mer_full_v,
        diff_u    = (bayes_v - mu_all_v) - (mer_full_v - mu_all_v),
        stringsAsFactors = FALSE
      )
    }))

    cat("\nFactor-level comparison by term:\n")
    cat(sprintf("  %s_raw   = coef(%s)\n", mer_label, mer_label))
    cat(sprintf(
      "  %s_full  = mu_all + (%s_raw - anchor)  [full b_j in lmebayes notation]\n",
      mer_label, mer_label
    ))
    cat(sprintf(
      "  u_%s     = %s_full - mu_all  (RE deviation from prior mean)\n",
      mer_label, mer_label
    ))
    cat(sprintf(
      "  u_%s     = %s - mu_all      (same deviation from posterior mean)\n",
      bayes_label, bayes_label
    ))
    cat(sprintf(
      "  diff_bf  = %s - %s_full;  diff_u = u_%s - u_%s (= diff_bf)\n\n",
      bayes_label, mer_label, bayes_label, mer_label
    ))

    num_long <- c(
      "mer_raw", "mer_full", "mu_all", "bayes",
      "u_mer", "u_bayes", "diff_bf", "diff_u"
    )
    out_long <- level_long
    out_long[num_long] <- lapply(out_long[num_long], function(x) round(x, digits))
    print(out_long)

    cat(sprintf(
      "\nMean |%s - %s_full| across factor levels, by term:\n",
      bayes_label, mer_label
    ))
    print(round(tapply(abs(level_long$diff_bf), level_long$term, mean), digits))
    cat(sprintf(
      "\nMean |u_%s - u_%s| (= |diff_bf|) across factor levels, by term:\n",
      bayes_label, mer_label
    ))
    print(round(tapply(abs(level_long$diff_u), level_long$term, mean), digits))
    cat(sprintf(
      "\nMean |%s - %s_raw| across factor levels, by term:\n",
      bayes_label, mer_label
    ))
    print(round(tapply(
      abs(level_long$bayes - level_long$mer_raw),
      level_long$term,
      mean
    ), digits))

    for (nm in re_names) {
      sub <- level_long[level_long$term == nm, , drop = FALSE]
      if (nrow(sub) < 2L) next
      cat(sprintf("\n--- %s ---\n", nm))
      if (nm %in% names(anchor)) {
        mu_ref <- if (nm %in% rownames(mu_mat)) {
          mu_mat[nm, 1L]
        } else {
          NA_real_
        }
        cat(sprintf(
          "  anchor = fixef %s: %s; mu_all (level 1): %s\n",
          nm,
          if (is.na(anchor[nm])) "NA" else sprintf("%.*f", digits, anchor[nm]),
          if (is.na(mu_ref)) "NA" else sprintf("%.*f", digits, mu_ref)
        ))
      }
      cat(sprintf(
        "  Cor(%s_full, %s): %s\n",
        mer_label, bayes_label,
        round(cor(sub$mer_full, sub$bayes), 3)
      ))
      cat(sprintf(
        "  Cor(u_%s, u_%s): %s\n",
        mer_label, bayes_label,
        round(cor(sub$u_mer, sub$u_bayes), 3)
      ))
      cat(sprintf(
        "  Mean |diff_bf|: %s\n",
        round(mean(abs(sub$diff_bf)), digits)
      ))
    }
  }

  invisible(list(
    mer_full = mer_full,
    bayes    = bayes,
    anchor   = anchor
  ))
}

#' Block~2 row kind for prior-aware \code{lmer} comparison
#' @noRd
.block2_fixef_row_kind <- function(component, parameter, prior_mean, mer_est) {
  if (identical(component, "(Intercept)") && identical(parameter, "(Intercept)")) {
    return("global_intercept")
  }
  if (is.finite(prior_mean) && abs(prior_mean) <= 1e-8) {
    return("null_effects")
  }
  if (is.finite(mer_est) && is.finite(prior_mean)) {
    tol <- max(1e-6, 0.05 * abs(mer_est))
    if (abs(prior_mean - mer_est) <= tol) {
      return("full_model_prior")
    }
  }
  "other"
}

#' Validate Block~2 fixed effects: MCMC vs ICM and prior-aware \code{lmer} checks
#'
#' Catches \code{simulate = FALSE} / mode-only fits (zero draw SD).  Under default
#' \code{Prior_Setup_GLMM(pwt = 0.01)} (\code{null_model} intercept,
#' \code{null_effects} slopes), non-intercept hyperparameters are expected to
#' shrink toward 0 relative to full-model \code{lmer} MLE; posterior SD is often below
#' \code{lmer} SE.  The global intercept prior mean comes from the null model
#' (reference in the table only); under weak \code{pwt} the posterior tracks
#' full-model \code{lmer}, not that prior mean.  Primary simulation check:
#' \code{|z_icm|} vs ICM / posterior mean.
#'
#' @param fit An \code{lmerb} or \code{glmerb} object with stored draws.
#' @param label Section label for messages.
#' @param z_icm_max Max \code{|z|} for \code{(mean - ICM) / SE(mean)}.
#' @param z_lmer_max Max \code{|z|} for \code{(mean - lmer) / lmer SE} on
#'   \code{full_model_prior} rows only.
#' @param shrink_slack Multiplier on \code{lmer SE} allowed beyond \code{|lmer est|}
#'   when checking null-effects attenuation toward 0.
#' @param se_ratio_min Lower bound on \code{post_sd / lmer SE} for
#'   \code{global_intercept} and \code{full_model_prior} rows (weak \code{pwt}:
#'   posterior SD should be near but typically below \code{lmer} SE).
#' @param se_ratio_max Upper bound on \code{post_sd / lmer SE} for all mapped rows.
#' @param se_ratio_min_null_effects Lower bound on \code{se_ratio} for
#'   \code{null_effects} rows. Defaults to \code{se_ratio_min}; set \code{NA} to skip
#'   (ING shrinkage may reduce posterior SD further).
#' @param check_se_ratio If \code{FALSE}, print \code{se_ratio} but do not stop
#'   (e.g. BlockEnvelopeCentering routes where \code{post_sd} is not comparable to
#'   \code{lmer} SE even when chain means track \code{lmer}).
#' @param digits Decimal places for printed table.
#' @return Invisibly, a data.frame of comparisons.
#' @keywords internal
#' @noRd
.validate_lmerb_block2_fixef_lmer <- function(
    fit,
    label = "lmerb",
    z_icm_max = 4,
    z_lmer_max = 4,
    shrink_slack = 0.5,
    se_ratio_min = 0.85,
    se_ratio_max = 1.05,
    se_ratio_min_null_effects = 0.85,
    check_se_ratio = TRUE,
    digits = 4L
) {
  if (!inherits(fit, c("lmerb", "glmerb"))) {
    stop("fit must be an lmerb or glmerb object.", call. = FALSE)
  }
  if (is.null(.lmerb_groupef_draws(fit))) {
    stop(
      "[", label, "] Block~2 validation requires simulate = TRUE ",
      "(fit$groupef is NULL).",
      call. = FALSE
    )
  }

  is_lmerb   <- inherits(fit, "lmerb")
  mer_label  <- if (is_lmerb) "lmer" else "glmer"
  mode_label <- if (is_lmerb) "ICM mean" else "post.mode"
  re_names    <- fit$model_setup$groupef.names
  popef_draws <- .lmerb_popef_draws(fit)
  popef_means <- .lmerb_popef_means(fit)
  popef_mode  <- .lmerb_popef_mode(fit)
  n_draws     <- nrow(popef_draws[[re_names[1L]]])
  mer_fit     <- .lmerb_reference_fit(fit)
  re_mod      <- fit$model_setup$popef.moderation

  rows <- lapply(re_names, function(k) {
    dm_k  <- popef_means[[k]]
    sd_k  <- apply(popef_draws[[k]], 2L, stats::sd)
    se_k  <- sd_k / sqrt(n_draws)
    icm_k <- popef_mode[[k]]
    pl_k <- fit$prior$pop.prior_list[[k]]
    lapply(names(dm_k), function(nm) {
      mer_ref <- .lmerb_lmer_fixef_lookup(
        mer_fit, k, nm, re_slope_moderation = re_mod
      )
      prior_mean <- unname(pl_k$mu[nm])
      data.frame(
        component   = k,
        parameter   = nm,
        prior_mean  = prior_mean,
        post_mean   = unname(dm_k[[nm]]),
        post_sd     = unname(sd_k[[nm]]),
        mc_se       = unname(se_k[[nm]]),
        icm         = unname(icm_k[[nm]]),
        mer_est     = mer_ref$estimate,
        mer_se      = mer_ref$se,
        stringsAsFactors = FALSE
      )
    })
  })
  cmp <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(cmp) <- NULL

  cmp$row_kind <- mapply(
    .block2_fixef_row_kind,
    cmp$component,
    cmp$parameter,
    cmp$prior_mean,
    cmp$mer_est,
    USE.NAMES = FALSE
  )
  cmp$z_icm <- with(cmp, (post_mean - icm) / mc_se)
  cmp$z_prior <- with(cmp, (post_mean - prior_mean) / mc_se)
  cmp$z_mer <- with(cmp, (post_mean - mer_est) / mer_se)
  cmp$se_ratio <- with(cmp, post_sd / mer_se)

  w_num <- digits + 4L
  fmt_num <- function(x) {
    formatC(x, digits = digits, width = w_num, format = "f", flag = " ")
  }

  cat(sprintf(
    "\n=== [%s] Block~2 fixed effects (n = %s draws) ===\n",
    label, format(n_draws, big.mark = ",")
  ))
  cat(
    "  Primary check: chain mean vs ICM (simulation; z_icm uses MC SE of the mean).  ",
    mer_label, " MLE is reference only;\n",
    if (check_se_ratio) {
      paste0(
        "  Uncertainty: post_sd vs ", mer_label,
        "_se via se_ratio (weak pwt: expect slightly below 1; enforced on mapped rows).\n"
      )
    } else {
      paste0(
        "  se_ratio (post_sd / ", mer_label,
        "_se) printed for reference only \u2014 not enforced on this route.\n"
      )
    },
    "  null_effects rows should shrink toward prior mean 0; ",
    "global intercept tracks full-model ", mer_label, " (null-model prior mean is reference only).\n\n",
    sep = ""
  )
  cat(sprintf(
    "  %-18s  %-22s  %*s  %*s  %*s  %*s  %*s  %*s  %*s  %6s  %6s\n",
    "RE component", "parameter",
    digits + 4L, "prior_mean",
    digits + 4L, paste0(mer_label, "_est"),
    digits + 4L, paste0(mer_label, "_se"),
    digits + 4L, "post_mean",
    digits + 4L, "post_sd",
    digits + 4L, "se_ratio",
    digits + 4L, mode_label,
    "z_icm", "z_mer"
  ))
  for (i in seq_len(nrow(cmp))) {
    r <- cmp[i, ]
    cat(
      sprintf("  %-18s  %-22s", r$component, r$parameter),
      fmt_num(r$prior_mean),
      fmt_num(r$mer_est),
      fmt_num(r$mer_se),
      fmt_num(r$post_mean),
      fmt_num(r$post_sd),
      if (is.finite(r$se_ratio)) fmt_num(r$se_ratio) else sprintf("%*s", w_num, "NA"),
      fmt_num(r$icm),
      sprintf("%6.2f", r$z_icm),
      if (is.finite(r$z_mer)) sprintf("%6.2f", r$z_mer) else "   NA",
      sprintf("  [%s]\n", r$row_kind),
      sep = "  "
    )
  }

  stopifnot(all(is.finite(cmp$post_mean)), all(is.finite(cmp$mc_se)))
  stopifnot(all(cmp$post_sd > 0), n_draws > 1L)

  mapped <- is.finite(cmp$mer_est) & is.finite(cmp$mer_se) & cmp$mer_se > 0

  .block2_stop_se_ratio <- function(bad, bound, threshold) {
    if (nrow(bad) == 0L) {
      return(invisible(NULL))
    }
    stop(
      "[", label, "] post_sd / ", mer_label, " SE ", bound, " ", threshold, " for: ",
      paste(
        paste0(
          bad$component, "::", bad$parameter,
          " (ratio=", sprintf("%.3f", bad$se_ratio), ")"
        ),
        collapse = "; "
      ),
      call. = FALSE
    )
  }

  se_mapped <- mapped & is.finite(cmp$se_ratio)
  if (isTRUE(check_se_ratio) && any(se_mapped)) {
    high <- cmp[se_mapped & cmp$se_ratio > se_ratio_max, , drop = FALSE]
    .block2_stop_se_ratio(high, ">", se_ratio_max)
    se_min_rows <- cmp$row_kind %in% c("global_intercept", "full_model_prior") &
      se_mapped
    low <- cmp[se_min_rows & cmp$se_ratio < se_ratio_min, , drop = FALSE]
    .block2_stop_se_ratio(low, "<", se_ratio_min)
    if (is.finite(se_ratio_min_null_effects)) {
      null_se <- cmp$row_kind == "null_effects" & se_mapped &
        cmp$se_ratio < se_ratio_min_null_effects
      low_null <- cmp[null_se, , drop = FALSE]
      .block2_stop_se_ratio(low_null, "<", se_ratio_min_null_effects)
    }
  }

  ## null_effects: attenuated toward 0 vs full-model lmer (same sign, |post| <= |lmer|)
  null_rows <- cmp$row_kind == "null_effects" & mapped
  if (any(null_rows)) {
    sub <- cmp[null_rows, , drop = FALSE]
    mer_sig <- abs(sub$mer_est) > sub$mer_se * 0.05
    if (any(mer_sig)) {
      sign_bad <- mer_sig & sign(sub$post_mean) != sign(sub$mer_est)
      if (any(sign_bad)) {
        bad <- sub[sign_bad, , drop = FALSE]
        stop(
          "[", label, "] null_effects sign mismatch vs ", mer_label, " for: ",
          paste(
            paste0(
              bad$component, "::", bad$parameter,
              " (post=", sprintf("%.4g", bad$post_mean),
              ", ", mer_label, "=", sprintf("%.4g", bad$mer_est), ")"
            ),
            collapse = "; "
          ),
          call. = FALSE
        )
      }
      shrink_bad <- mer_sig &
        (abs(sub$post_mean) > abs(sub$mer_est) + shrink_slack * sub$mer_se)
      if (any(shrink_bad)) {
        bad <- sub[shrink_bad, , drop = FALSE]
        stop(
          "[", label, "] null_effects not attenuated toward 0 (|post| > |",
          mer_label, "|) for: ",
          paste(
            paste0(
              bad$component, "::", bad$parameter,
              " (post=", sprintf("%.4g", bad$post_mean),
              ", ", mer_label, "=", sprintf("%.4g", bad$mer_est), ")"
            ),
            collapse = "; "
          ),
          call. = FALSE
        )
      }
    }
  }

  ## global_intercept: weak pwt => posterior near full-model lmer (not null prior mean)
  gi_rows <- cmp$row_kind == "global_intercept" & mapped
  if (any(gi_rows)) {
    if (any(abs(cmp$z_mer[gi_rows]) >= z_lmer_max, na.rm = TRUE)) {
      bad <- cmp[gi_rows & abs(cmp$z_mer) >= z_lmer_max, , drop = FALSE]
      stop(
        "[", label, "] |z_mer| >= ", z_lmer_max, " on global intercept for: ",
        paste(
          paste0(
            bad$component, "::", bad$parameter,
            " (z=", sprintf("%.2f", bad$z_mer), ")"
          ),
          collapse = "; "
        ),
        call. = FALSE
      )
    }
  }

  ## full_model_prior / other: symmetric z vs lmer MLE
  lmer_sym <- cmp$row_kind %in% c("full_model_prior", "other") & mapped
  if (any(lmer_sym)) {
    if (any(abs(cmp$z_mer[lmer_sym]) >= z_lmer_max, na.rm = TRUE)) {
      bad <- cmp[lmer_sym & abs(cmp$z_mer) >= z_lmer_max, , drop = FALSE]
      stop(
        "[", label, "] |z_mer| >= ", z_lmer_max, " vs ", mer_label, " for: ",
        paste(
          paste0(bad$component, "::", bad$parameter, " (z=", sprintf("%.2f", bad$z_mer), ")"),
          collapse = "; "
        ),
        call. = FALSE
      )
    }
  }

  if (any(abs(cmp$z_icm) >= z_icm_max, na.rm = TRUE)) {
    bad <- cmp[which(abs(cmp$z_icm) >= z_icm_max), , drop = FALSE]
    stop(
      "[", label, "] |z_icm| >= ", z_icm_max, " for: ",
      paste(
        paste0(bad$component, "::", bad$parameter, " (z=", sprintf("%.2f", bad$z_icm), ")"),
        collapse = "; "
      ),
      call. = FALSE
    )
  }

  invisible(cmp)
}

#' Manual-test defaults: \code{post_sd / mer_se} on all mapped Block~2 rows.
#' @noRd
BLOCK2_FIXEF_SE <- list(
  se_ratio_min = 0.85,
  se_ratio_max = 1.05
)

#' ING manual tests: no lower bound on \code{null_effects} \code{se_ratio}.
#' @noRd
BLOCK2_FIXEF_SE_ING <- list(
  se_ratio_min = 0.85,
  se_ratio_max = 1.05,
  se_ratio_min_null_effects = NA_real_
)

#' Block~2 fixef validation wrapper for manual MER scripts (\code{pwt = 0.01}).
#' @noRd
.validate_manual_block2_fixef <- function(
    fit,
    label,
    z_icm_max = 4,
    ing = FALSE,
    check_se_ratio = TRUE,
    ...
) {
  se_args <- if (isTRUE(ing)) BLOCK2_FIXEF_SE_ING else BLOCK2_FIXEF_SE
  do.call(
    .validate_lmerb_block2_fixef_lmer,
    c(
      list(
        fit = fit,
        label = label,
        z_icm_max = z_icm_max,
        check_se_ratio = check_se_ratio
      ),
      se_args,
      list(...)
    )
  )
}

#' Validate lmerb random effects: print cor / table, then ICM and lmer_full thresholds
#' @noRd
.validate_lmerb_re <- function(
    fit,
    label = "lmerb",
    cor_icm_min = 0.9,
    cor_full_min = 0.85,
    mode_match_min = 0.9
) {
  stopifnot(inherits(fit, "lmerb"))
  re_names <- fit$model_setup$groupef.names
  grp_col  <- fit$model_setup$group_name
  grp_levs <- rownames(
    lmebayesCore:::.lmebayes_reference_coef(.lmerb_reference_fit(fit))[[grp_col]]
  )
  J        <- length(grp_levs)
  icm_b    <- .lmerb_groupef_mode(fit)
  groupef  <- .lmerb_groupef_draws(fit)
  n_draws  <- nrow(.lmerb_popef_draws(fit)[[re_names[1L]]])

  stopifnot(setequal(rownames(icm_b), grp_levs))
  stopifnot(setequal(colnames(.lmerb_popef_mu(fit)), grp_levs))
  stopifnot(setequal(unique(as.character(groupef[[grp_col]])), grp_levs))
  stopifnot(identical(
    sort(table(as.character(groupef[[grp_col]]))),
    sort(table(rep(grp_levs, n_draws)))
  ))
  cat(sprintf("\n[%s] Ordering OK across groupef.mode / mu_all / groupef\n", label))

  re_draws_mean <- tapply(
    seq_len(nrow(groupef)),
    groupef[[grp_col]],
    function(idx) colMeans(groupef[idx, re_names, drop = FALSE]),
    simplify = FALSE
  )

  mcmc_mat <- do.call(rbind, lapply(grp_levs, function(l) {
    re_draws_mean[[as.character(l)]]
  }))
  rownames(mcmc_mat) <- grp_levs
  icm_mat <- icm_b[grp_levs, re_names, drop = FALSE]

  n_match <- NA_integer_
  if (length(re_names) > 1L) {
    scl <- apply(icm_mat, 2L, sd)
    scl[scl < 1e-8] <- 1
    A <- sweep(mcmc_mat, 2L, scl, "/")
    B <- sweep(icm_mat, 2L, scl, "/")
    D <- outer(rowSums(A^2), rep(1, J)) +
      outer(rep(1, J), rowSums(B^2)) - 2 * A %*% t(B)
    nearest <- apply(D, 1L, which.min)
    n_match <- sum(nearest == seq_len(J))
    cat(sprintf(
      "[%s] Nearest-mode matching: %d of %d groups on own ICM mode\n",
      label, n_match, J
    ))
  } else {
    cat(sprintf(
      "[%s] Nearest-mode matching: skipped (univariate RE; use cor vs own ICM mode)\n",
      label
    ))
  }

  mer_full <- .mer_re_reference_full(fit)
  cat(sprintf("[%s] cor(MCMC mean vs ICM / lmer_full):\n", label))
  for (k in re_names) {
    c_icm  <- cor(mcmc_mat[, k], icm_mat[, k])
    c_full <- cor(mcmc_mat[, k], mer_full[, k])
    cat(sprintf("  %-14s  vs ICM: %6.3f   vs lmer_full: %6.3f\n", k, c_icm, c_full))
  }

  cat(sprintf("\n=== Random effects table (%s) ===\n\n", label))
  print_mer_bayes_re_compare(fit)

  if (length(re_names) > 1L) {
    stopifnot(n_match >= ceiling(mode_match_min * J))
  }
  for (k in re_names) {
    c_icm  <- cor(mcmc_mat[, k], icm_mat[, k])
    c_full <- cor(mcmc_mat[, k], mer_full[, k])
    stopifnot(c_icm >= cor_icm_min)
    stopifnot(c_full >= cor_full_min)
  }

  invisible(list(mcmc = mcmc_mat, mer_full = mer_full))
}

#' Validate glmerb random effects: print cor / table, then ICM and glmer_full thresholds
#' @noRd
.validate_glmerb_re <- function(
    fit,
    label = "glmerb",
    cor_icm_min = 0.9,
    cor_full_min = 0.8,
    mode_match_min = 0.9
) {
  stopifnot(inherits(fit, "glmerb"))
  re_names <- fit$model_setup$groupef.names
  grp_col  <- fit$model_setup$group_name
  grp_levs <- rownames(
    lmebayesCore:::.lmebayes_reference_coef(.lmerb_reference_fit(fit))[[grp_col]]
  )
  J        <- length(grp_levs)
  icm_b    <- .lmerb_groupef_mode(fit)
  groupef  <- .lmerb_groupef_draws(fit)
  n_draws  <- nrow(.lmerb_popef_draws(fit)[[re_names[1L]]])

  stopifnot(setequal(rownames(icm_b), grp_levs))
  stopifnot(setequal(colnames(.lmerb_popef_mu(fit)), grp_levs))
  stopifnot(setequal(unique(as.character(groupef[[grp_col]])), grp_levs))
  stopifnot(identical(
    sort(table(as.character(groupef[[grp_col]]))),
    sort(table(rep(grp_levs, n_draws)))
  ))
  cat(sprintf("\n[%s] Ordering OK across groupef.mode / mu_all / groupef\n", label))

  re_draws_mean <- tapply(
    seq_len(nrow(groupef)),
    groupef[[grp_col]],
    function(idx) colMeans(groupef[idx, re_names, drop = FALSE]),
    simplify = FALSE
  )

  mcmc_mat <- do.call(rbind, lapply(grp_levs, function(l) {
    re_draws_mean[[as.character(l)]]
  }))
  rownames(mcmc_mat) <- grp_levs
  icm_mat <- icm_b[grp_levs, re_names, drop = FALSE]

  n_match <- NA_integer_
  if (length(re_names) > 1L) {
    scl <- apply(icm_mat, 2L, sd)
    scl[scl < 1e-8] <- 1
    A <- sweep(mcmc_mat, 2L, scl, "/")
    B <- sweep(icm_mat, 2L, scl, "/")
    D <- outer(rowSums(A^2), rep(1, J)) +
      outer(rep(1, J), rowSums(B^2)) - 2 * A %*% t(B)
    nearest <- apply(D, 1L, which.min)
    n_match <- sum(nearest == seq_len(J))
    cat(sprintf(
      "[%s] Nearest-mode matching: %d of %d groups on own ICM mode\n",
      label, n_match, J
    ))
  } else {
    cat(sprintf(
      "[%s] Nearest-mode matching: skipped (univariate RE; use cor vs own ICM mode)\n",
      label
    ))
  }

  mer_full <- .mer_re_reference_full(fit)
  cat(sprintf("[%s] cor(MCMC mean vs ICM / glmer_full):\n", label))
  for (k in re_names) {
    c_icm  <- cor(mcmc_mat[, k], icm_mat[, k])
    c_full <- cor(mcmc_mat[, k], mer_full[, k])
    cat(sprintf("  %-14s  vs ICM: %6.3f   vs glmer_full: %6.3f\n", k, c_icm, c_full))
  }

  cat(sprintf("\n=== Random effects table (%s) ===\n\n", label))
  print_mer_bayes_re_compare(fit)

  if (length(re_names) > 1L) {
    stopifnot(n_match >= ceiling(mode_match_min * J))
  }
  for (k in re_names) {
    c_icm  <- cor(mcmc_mat[, k], icm_mat[, k])
    c_full <- cor(mcmc_mat[, k], mer_full[, k])
    stopifnot(c_icm >= cor_icm_min)
    stopifnot(c_full >= cor_full_min)
  }

  invisible(list(mcmc = mcmc_mat, mer_full = mer_full))
}

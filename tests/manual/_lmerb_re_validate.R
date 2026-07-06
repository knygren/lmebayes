#' Shared lmerb RE validation (ordering + ICM / lmer_full cor + print table).
#'
#' Primary check: MCMC mean vs \code{ranef.mode} (simulation).  \code{lmer_full}
#' is a reference alignment only; under default priors cor can sit below 1 even
#' when the sampler is healthy (MC noise at moderate \code{n}).
#'
#' @param fit A fitted \code{lmerb} object with stored \code{coefficients}.
#' @param label Section label for messages.
#' @param cor_icm_min Minimum cor(MCMC mean, ICM mode) per RE term.
#' @param cor_full_min Minimum cor(MCMC mean, lmer_full) per RE term.
#' @param mode_match_min Minimum fraction of groups whose chain mean is nearest
#'   its own ICM mode vector.  Only applied when \code{p_re > 1} (multivariate
#'   RE); skipped for intercept-only models where scalar modes can collide.
#' @keywords internal
.validate_lmerb_re <- function(
    fit,
    label = "lmerb",
    cor_icm_min = 0.9,
    cor_full_min = 0.85,
    mode_match_min = 0.9
) {
  stopifnot(inherits(fit, "lmerb"))
  re_names <- fit$model_setup$re_coef_names
  grp_col  <- fit$model_setup$group_name
  grp_levs <- rownames(coef(fit$lmer)[[grp_col]])
  J        <- length(grp_levs)
  icm_b    <- fit$ranef.mode
  n_draws  <- nrow(fit$fixef[[re_names[1L]]])

  stopifnot(setequal(rownames(icm_b), grp_levs))
  stopifnot(setequal(colnames(fit$fixef.mu), grp_levs))
  stopifnot(setequal(unique(as.character(fit$coefficients[[grp_col]])), grp_levs))
  stopifnot(identical(
    sort(table(as.character(fit$coefficients[[grp_col]]))),
    sort(table(rep(grp_levs, n_draws)))
  ))
  cat(sprintf("\n[%s] Ordering OK across ranef.mode / mu_all / coefficients\n", label))

  re_draws_mean <- tapply(
    seq_len(nrow(fit$coefficients)),
    fit$coefficients[[grp_col]],
    function(idx) colMeans(fit$coefficients[idx, re_names, drop = FALSE]),
    simplify = FALSE
  )

  mcmc_mat <- do.call(rbind, lapply(grp_levs, function(l) {
    re_draws_mean[[as.character(l)]]
  }))
  rownames(mcmc_mat) <- grp_levs
  icm_mat <- icm_b[grp_levs, re_names, drop = FALSE]

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
    stopifnot(n_match >= ceiling(mode_match_min * J))
  } else {
    cat(sprintf(
      "[%s] Nearest-mode matching: skipped (univariate RE; use cor vs own ICM mode)\n",
      label
    ))
  }

  mer_full <- lmebayes:::.mer_re_reference_full(fit)
  cat(sprintf("[%s] cor(MCMC mean vs ICM / lmer_full):\n", label))
  for (k in re_names) {
    c_icm  <- cor(mcmc_mat[, k], icm_mat[, k])
    c_full <- cor(mcmc_mat[, k], mer_full[, k])
    cat(sprintf("  %-14s  vs ICM: %6.3f   vs lmer_full: %6.3f\n", k, c_icm, c_full))
    stopifnot(c_icm >= cor_icm_min)
    stopifnot(c_full >= cor_full_min)
  }

  cat(sprintf("\n=== Random effects table (%s) ===\n\n", label))
  lmebayes:::print_mer_bayes_re_compare(fit)
  invisible(list(mcmc = mcmc_mat, mer_full = mer_full))
}

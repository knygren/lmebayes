# Smoke test: glmerb (Poisson) on bayesrules::airbnb with level-2 covariates
# and cross-level RE moderation (walk_c:rating_c, transit_c:log_price_c).
#
#   Rscript data-raw/test_glmerb_airbnb.R
#   Rscript data-raw/test_glmerb_airbnb.R quick
#   Rscript data-raw/test_glmerb_airbnb.R small

args <- commandArgs(trailingOnly = TRUE)
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))
use_small <- any(tolower(args) %in% c("small", "--small"))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

if (use_small) {
  data("airbnb_small", package = "bayesrules")
  dat <- airbnb_small
  message("Using bayesrules::airbnb_small (n = ", nrow(dat), ")")
} else {
  data("airbnb", package = "bayesrules")
  dat <- airbnb
  message("Using bayesrules::airbnb (n = ", nrow(dat), ")")
}

dat$rating_c    <- dat$rating - mean(dat$rating)
dat$log_price_c <- scale(log(dat$price + 1))[, 1]
dat$walk_c      <- dat$walk_score - mean(dat$walk_score)
dat$transit_c   <- dat$transit_score - mean(dat$transit_score)
dat <- dat[complete.cases(dat[, c(
  "reviews", "rating", "rating_c", "price",
  "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
)]), ]

form <- reviews ~
  walk_c + transit_c +
  rating_c + log_price_c +
  walk_c:rating_c + transit_c:log_price_c +
  (1 + rating_c + log_price_c || neighborhood)

design <- model_setup(form, data = dat, family = poisson())
stopifnot(isTRUE(design$rank_ok))
stopifnot(length(design$re_coef_names) == 3L)
stopifnot(ncol(design$X_hyper[["(Intercept)"]]) >= 3L)
stopifnot(ncol(design$X_hyper[["rating_c"]]) >= 2L)
stopifnot(ncol(design$X_hyper[["log_price_c"]]) >= 2L)

ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)
stopifnot(inherits(ps, "lmebayes_prior_setup"))
stopifnot(is.null(ps$dispersion_ranef))
stopifnot(all(diag(ps$Sigma_ranef) > 0))

n_draw <- if (run_quick) 50L else 200L
message("Posterior draws: n = ", n_draw)

fit <- glmerb(
  form,
  data = dat,
  family = poisson(),
  pfamily_list = pfamily_list(ps),
  n = n_draw,
  seed = 42L
)

stopifnot(inherits(fit, "glmerb"))
stopifnot(length(fit$fixef.mode) == 3L)
re_names <- fit$model_setup$re_coef_names
stopifnot(nrow(fit$fixef[[re_names[1L]]]) == n_draw)

print(summary(fit))

## ===========================================================================
## Ordering diagnostics (Poisson / GLM path)
##
## Mirror of the per-group z-test in test_glmerb_big_word_club.R, with one
## caveat: for a Poisson model the posterior is skewed, so the MCMC mean can
## differ from the ICM *mode* even when the sampler is correct.  The z-table
## below is therefore informational.  The hard assertions are the ordering
## checks: a scrambled neighborhood ordering would make each group's MCMC
## mean land near some *other* group's ICM mode, which the nearest-mode
## matching test and the cross-group correlations detect immediately.
## ===========================================================================

grp_col  <- fit$model_setup$group_name
grp_levs <- rownames(coef(fit$glmer)[[grp_col]])
J        <- length(grp_levs)
icm_b    <- fit$ranef.mode   # J x p_re, rownames = group levels
n_draws  <- nrow(fit$fixef[[re_names[1L]]])

## --- Label consistency across all returned objects -------------------------
stopifnot(setequal(rownames(icm_b), grp_levs))
stopifnot(setequal(colnames(fit$fixef.mu), grp_levs))
stopifnot(setequal(unique(as.character(fit$coefficients[[grp_col]])), grp_levs))
stopifnot(identical(sort(table(as.character(fit$coefficients[[grp_col]]))),
                    sort(table(rep(grp_levs, n_draws)))))
cat("\nOrdering: group labels consistent across ranef.mode / mu_all / coefficients\n")

## --- Per-group MCMC means and SDs (indexed by level name, not position) ----
re_draws_mean <- tapply(
  seq_len(nrow(fit$coefficients)),
  fit$coefficients[[grp_col]],
  function(idx) colMeans(fit$coefficients[idx, re_names, drop = FALSE]),
  simplify = FALSE
)
re_draws_sd <- tapply(
  seq_len(nrow(fit$coefficients)),
  fit$coefficients[[grp_col]],
  function(idx) apply(fit$coefficients[idx, re_names, drop = FALSE], 2L, sd),
  simplify = FALSE
)

## --- Informational z-table: MCMC mean vs ICM mode --------------------------
cat("\n=== Random effects: MCMC mean vs ICM mode (informational for Poisson) ===\n\n")
cat(sprintf("  %-22s  %-14s  %10s  %10s  %10s  %6s\n",
            "group", "RE component", "MCMC mean", "ICM mode", "SE(mean)", "z"))
cat(sprintf("  %-22s  %-14s  %10s  %10s  %10s  %6s\n",
            strrep("-", 22L), strrep("-", 14L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L), strrep("-", 6L)))
n_flagged <- 0L
for (lev in grp_levs) {
  lev_chr <- as.character(lev)
  for (k in re_names) {
    mcmc_m <- re_draws_mean[[lev_chr]][[k]]
    mcmc_s <- re_draws_sd[[lev_chr]][[k]]
    icm_m  <- icm_b[lev_chr, k]
    se_val <- mcmc_s / sqrt(n_draws)
    z_val  <- (mcmc_m - icm_m) / se_val
    flag   <- if (abs(z_val) > 3) " *" else "  "
    if (abs(z_val) > 3) n_flagged <- n_flagged + 1L
    cat(sprintf("  %-22s  %-14s  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                lev_chr, k, mcmc_m, icm_m, se_val, z_val, flag))
  }
}
total_tests <- J * length(re_names)
cat(sprintf(
  "\n  %d of %d tests flagged |z| > 3  (mean-vs-mode gaps are expected under Poisson skew)\n",
  n_flagged, total_tests
))

## --- Hard check 1: nearest-mode matching (scramble detector) ---------------
## Each group's MCMC mean vector must be nearest (standardized Euclidean) to
## its OWN ICM mode.  Skewness shifts every group by a small, similar amount;
## a permuted ordering moves groups onto other groups' modes.
mcmc_mat <- do.call(rbind, lapply(grp_levs, function(l) re_draws_mean[[as.character(l)]]))
rownames(mcmc_mat) <- grp_levs
icm_mat <- icm_b[grp_levs, re_names, drop = FALSE]

scl <- apply(icm_mat, 2L, sd)
scl[scl < 1e-8] <- 1
A <- sweep(mcmc_mat, 2L, scl, "/")
B <- sweep(icm_mat, 2L, scl, "/")
D <- outer(rowSums(A^2), rep(1, J)) + outer(rep(1, J), rowSums(B^2)) - 2 * A %*% t(B)
nearest <- apply(D, 1L, which.min)
n_match <- sum(nearest == seq_len(J))
cat(sprintf(
  "\nNearest-mode matching: %d of %d groups matched to their own ICM mode\n",
  n_match, J
))
if (n_match < J) {
  mism <- which(nearest != seq_len(J))
  for (i in mism) {
    cat(sprintf("  MISMATCH: MCMC mean of '%s' nearest to ICM mode of '%s'\n",
                grp_levs[i], grp_levs[nearest[i]]))
  }
}
stopifnot(n_match >= ceiling(0.9 * J))

## --- Hard check 2: cross-group correlations --------------------------------
## Across neighborhoods, MCMC means must track both the ICM modes and the
## independent lme4 reference coef(glmer).  Scrambled ordering drives these
## toward zero.
glmer_mat <- as.matrix(coef(fit$glmer)[[grp_col]][grp_levs, re_names, drop = FALSE])
cat("\nCross-group correlations (MCMC mean vs ICM mode, MCMC mean vs coef(glmer)):\n")
for (k in re_names) {
  c_icm   <- cor(mcmc_mat[, k], icm_mat[, k])
  c_glmer <- cor(mcmc_mat[, k], glmer_mat[, k])
  cat(sprintf("  %-14s  vs ICM mode: %6.3f   vs coef(glmer): %6.3f\n",
              k, c_icm, c_glmer))
  stopifnot(c_icm > 0.9)
  stopifnot(c_glmer > 0.8)
}

cat("\ntest_glmerb_airbnb: OK\n")

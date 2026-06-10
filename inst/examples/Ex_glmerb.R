## glmerb: neighborhood random effects on bayesrules::airbnb
##
## Poisson model for review counts. Listing-level predictors (rating, log price)
## have fixed main effects and neighborhood random slopes. Neighborhood-level
## walkability and transit scores are level-2 covariates on the random intercept;
## they also moderate the random slopes via cross-level interactions
## (walk_c:rating_c, transit_c:log_price_c), analogous to
## free_reduced_lunch:distracted_a1 in inst/examples/Ex_lmerb.R.
##
## Workflow: model_setup(), then Prior_Setup_lmebayes(), then glmerb().

data(airbnb, package = "bayesrules")

dat <- airbnb
dat$rating_c    <- dat$rating - mean(dat$rating)
dat$log_price_c <- scale(log(dat$price + 1))[, 1]
dat$walk_c      <- dat$walk_score - mean(dat$walk_score)
dat$transit_c   <- dat$transit_score - mean(dat$transit_score)
dat <- dat[complete.cases(dat[, c(
  "reviews", "rating", "rating_c", "price", "log_price_c",
  "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
)]), ]

form_glmer <- reviews ~
  walk_c + transit_c +
  rating_c + log_price_c +
  walk_c:rating_c + transit_c:log_price_c +
  (1 + rating_c + log_price_c || neighborhood)

design <- model_setup(form_glmer, data = dat, family = poisson())
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(form_glmer, data = dat, family = poisson(), pwt = 0.01)
cat("\n=== Prior_Setup_lmebayes ===\n\n")
print(ps)

fit <- glmerb(
  form_glmer,
  data = dat,
  family = poisson(),
  measurement_prior_list = ps,
  n = 1000L,
  seed = 42L
)

cat("\n=== summary(fit) ===\n\n")
print(summary(fit))

re_names <- fit$model_setup$re_coef_names
n_draws  <- nrow(fit$fixef_draws[[re_names[1L]]])

## --- Block 2 fixed effects: MCMC mean vs ICM posterior mode ------------------
##
## Block 2 updates are multivariate-normal (Gaussian regression of b_k on
## X_hyper). The ICM vector coef.mode is the exact joint mode given fixed
## variance components; colMeans(fixef_draws) is the MCMC estimate of the
## posterior mean (which may differ from the mode in Poisson models).
##
## Per coefficient: z = (draws mean - ICM mode) / SE(mean), SE = SD/sqrt(n).
## Combined (stacked gamma): Hotelling T^2 on the p-dimensional draw mean
## vs coef.mode. Under i.i.d. MVN draws, T^2 maps to an F(p, n-p) (equivalently
## chi-sq with df = p for large n). This is the multivariate analogue of an
## F-test for regression coefficients; it accounts for correlation among the
## hyperparameters within a draw. A sum of z^2 (independence) is reported as
## a simpler chi-sq reference when p is moderate.

block2_draw_matrix <- function(fit) {
  re <- fit$model_setup$re_coef_names
  parts <- lapply(re, function(k) {
    m <- fit$fixef_draws[[k]]
    colnames(m) <- paste0(k, "::", colnames(m))
    m
  })
  do.call(cbind, parts)
}

block2_mode_vector <- function(fit) {
  re <- fit$model_setup$re_coef_names
  unlist(lapply(re, function(k) {
    v <- fit$coef.mode[[k]]
    names(v) <- paste0(k, "::", names(v))
    v
  }))
}

hotelling_test <- function(draw_mat, target_vec) {
  mcmc_mean <- colMeans(draw_mat)
  diff_vec  <- mcmc_mean - target_vec
  n         <- nrow(draw_mat)
  p         <- ncol(draw_mat)
  if (p >= n) {
    return(list(
      p = p, n = n, diff = diff_vec, T2 = NA_real_, F = NA_real_,
      p_F = NA_real_, p_chisq = NA_real_, p_z2 = NA_real_,
      note = "Need n > p for Hotelling F-test."
    ))
  }
  S <- stats::cov(draw_mat)
  Sinv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  T2 <- as.numeric(n * t(diff_vec) %*% Sinv %*% diff_vec)
  df2 <- n - p
  F_stat <- df2 / (p * max(n - 1L, 1L)) * T2
  se <- apply(draw_mat, 2L, sd) / sqrt(n)
  z <- diff_vec / se
  list(
    p = p, n = n, diff = diff_vec, z = z,
    T2 = T2,
    F = F_stat,
    p_F = stats::pf(F_stat, df1 = p, df2 = df2, lower.tail = FALSE),
    p_chisq = stats::pchisq(T2, df = p, lower.tail = FALSE),
    p_z2 = stats::pchisq(sum(z^2), df = p, lower.tail = FALSE),
    note = NULL
  )
}

cat("\n=== Block 2 fixed effects: MCMC mean vs ICM posterior mode ===\n\n")
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            "RE component", "parameter", "draws mean", "draws SD", "SE(mean)",
            "ICM mode", "z"))
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            strrep("-", 18L), strrep("-", 28L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L),
            strrep("-", 10L), strrep("-", 7L)))

n_flagged <- 0L
total_fe  <- 0L
for (k in re_names) {
  dm_k  <- fit$coef.means[[k]]
  sd_k  <- apply(fit$fixef_draws[[k]], 2L, sd)
  se_k  <- sd_k / sqrt(n_draws)
  icm_k <- fit$coef.mode[[k]]
  for (nm in names(dm_k)) {
    z_val <- (dm_k[[nm]] - icm_k[[nm]]) / se_k[[nm]]
    total_fe <- total_fe + 1L
    flag  <- if (abs(z_val) > 2) " *" else "  "
    if (abs(z_val) > 2) n_flagged <- n_flagged + 1L
    cat(sprintf("  %-18s  %-28s  %10.4f  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                k, nm, dm_k[[nm]], sd_k[[nm]], se_k[[nm]], icm_k[[nm]], z_val, flag))
  }
}
cat(sprintf(
  "\n  %d of %d coefficients flagged |z| > 2  (expected ~%.1f by chance at 5%% two-sided)\n",
  n_flagged, total_fe, total_fe * 0.05
))
cat("  (* |z| > 2: marginal MCMC mean differs from ICM mode)\n")

cat("\n=== Combined tests (Hotelling T^2  ~>  F; chi-sq references) ===\n\n")

draw_mat <- block2_draw_matrix(fit)
mode_vec <- block2_mode_vector(fit)
glob <- hotelling_test(draw_mat, mode_vec)

cat(sprintf(
  "  All Block 2 coefficients (p = %d, n = %d draws):\n",
  glob$p, glob$n
))
cat(sprintf("    || MCMC mean - ICM mode ||_2 = %.4f\n",
            sqrt(sum(glob$diff^2))))
if (is.na(glob$T2)) {
  cat("   ", glob$note, "\n")
} else {
  cat(sprintf("    Hotelling T^2 = %.4f\n", glob$T2))
  cat(sprintf("    F(%d, %d) = %.4f,  p = %.4g  (preferred combined test)\n",
              glob$p, glob$n - glob$p, glob$F, glob$p_F))
  cat(sprintf("    Chi-sq(df = %d) on T^2:        p = %.4g  (large-n reference)\n",
              glob$p, glob$p_chisq))
  cat(sprintf("    Sum of z^2 (independence):     p = %.4g  (ignores within-draw correlation)\n",
              glob$p_z2))
}

cat("\n  Per RE component (separate Hotelling on each gamma_k block):\n")
for (k in re_names) {
  dm_k <- fit$fixef_draws[[k]]
  icm_k <- fit$coef.mode[[k]]
  ht <- hotelling_test(dm_k, icm_k)
  cat(sprintf(
    "    [%s]  q = %d:  F = %.3f, p = %.4g  ||diff|| = %.4f\n",
    k, ht$p, ht$F, ht$p_F, sqrt(sum(ht$diff^2))
  ))
}

cat(paste0(
  "\n  Interpretation: these are MCMC-vs-ICM diagnostics, not classical ",
  "hypothesis tests.\n  Poisson Block 1 (non-Gaussian) can make the posterior ",
  "mean differ from the ICM\n  mode even when the sampler is well mixed. ",
  "Use n >= 500+ for stable combined p-values.\n"
))

grp_col <- fit$model_setup$group_name
cat("\n=== mu_all varies by neighborhood (rating_c slope prior mean) ===\n\n")
mu_rating <- fit$mu_all["rating_c", , drop = TRUE]
walk_by_nbhd <- tapply(dat$walk_c, dat$neighborhood, function(x) x[1L])
walk_by_nbhd <- walk_by_nbhd[names(mu_rating)]
cor_mu_walk <- cor(mu_rating, walk_by_nbhd, use = "complete.obs")
cat(sprintf(
  "  Cor(mu_all[rating_c], neighborhood walk_c): %.3f\n",
  cor_mu_walk
))
cat("  (Positive correlation: higher walkability => higher prior mean rating slope.)\n")

cat("\nNeighborhood random effects (first 6 levels):\n")
ri_glmer <- coef(fit$glmer)[[grp_col]][, re_names, drop = FALSE]
ri_mode  <- fit$ranef.mode[, re_names, drop = FALSE]
print(head(cbind(glmer = ri_glmer, glmerb = ri_mode)))

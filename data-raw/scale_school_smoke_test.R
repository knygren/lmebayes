## Scaling smoke test: Block~1 ING dGamma with increasing number of schools.
##
## Uses .prepare_small5_lmerb_manual(n_schools = k) which picks the k largest
## full-rank schools from big_word_club.  Runs n = 50 draws at each combination
## of k ∈ {10, 15, 20, 25} x pwt_measurement ∈ {0.49, 0.20, 0.01}.
## For each cell prints: window, timing, ranef.iters.mean, sigma2.mean.
## Aborts a cell on error but continues to the next.

pkg_root <- local({
  p <- normalizePath(getwd())
  while (!file.exists(file.path(p, "DESCRIPTION")) && dirname(p) != p) p <- dirname(p)
  p
})
setwd(pkg_root)
source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

## ---------------------------------------------------------------------------
## EXPLORATORY: implied per-school (group-level) dispersion estimates.
##
## dispersion_ranef (sigma^2) is currently one pooled REML estimate shared
## across all schools (design$residual_var = sigma(lmer_fit)^2). As a first
## step toward exploring school-specific sigma^2_j priors informed by a
## classical group-level estimate, compute a per-school dispersion estimate
## from the SAME reference lmer fit that Prior_Setup_lmebayes/model_setup
## already builds (all full-rank schools) and compare it to the pooled
## estimate. No new production function -- exploratory script code only.
## ---------------------------------------------------------------------------
fx_full <- .prepare_small5_all_full_rank_manual()

mer_fit  <- fx_full$design$lmer_fit        # exact reference fit model_setup() builds
re_names <- fx_full$design$re_coef_names   # e.g. c("(Intercept)", "distracted_ppvt")
grp_col  <- fx_full$design$group_name
groups   <- fx_full$design$groups          # per-observation factor, same order as residuals()
p_re     <- length(re_names)

stopifnot(identical(nlevels(groups), nrow(coef(mer_fit)[[grp_col]])))

## Per-school random-effects estimates: "full" per-group coefficients
## (fixef + BLUP), i.e. the same b_j representation lmerb's Block~1 draws
## target -- same columns as re_coef_names, consistent with how lmerb reports
## random effects (see .mer_re_reference_full() in
## R/lmebayes_mer_re_compare.R for the analogous lmerb-fit-based version).
b_j <- as.matrix(coef(mer_fit)[[grp_col]][, re_names, drop = FALSE])
b_j <- b_j[levels(groups), , drop = FALSE]

## Classical per-school dispersion: RSS_j / (n_j - p_re) from the response
## residuals of the SAME fit (fixef + BLUP fitted values), grouped by school.
resid_all <- stats::residuals(mer_fit)
n_j       <- as.integer(table(groups))
rss_j     <- as.numeric(tapply(resid_all, groups, function(r) sum(r^2)))
names(n_j)   <- levels(groups)
names(rss_j) <- levels(groups)

sigma2_hat_j <- ifelse(n_j > p_re, rss_j / (n_j - p_re), NA_real_)
sigma2_pooled <- fx_full$design$residual_var   # = sigma(mer_fit)^2

## Per-school ING *window* shape/rate: the limiting-posterior Gamma on
## 1/sigma2_j mean-matched to THIS school's own classical sigma2_hat_j, with
## shape set by this school's own df_j -- same convention
## ING_TRUNCATION_WINDOW.md already uses for tau2_k (a_inf = (J+1)/2,
## b_inf = tau2*(J-1)/2), applied here per school instead of per RE
## component. This is NOT the dGamma() prior actually fed to the sampler
## (see inst/GROUP_DISPERSION_HYPERPRIOR.md Sec.8 -- that stays the SHARED
## (a0, b0) hyperprior computed below); it is only a per-school reference
## for how wide/where that school's own envelope truncation window would be
## centered if computed independently of the pooled fit.
df_j_all      <- n_j - p_re
prior_shape_j <- ifelse(df_j_all > 1, (df_j_all + 1) / 2, NA_real_)
prior_rate_j  <- ifelse(df_j_all > 1, sigma2_hat_j * (df_j_all - 1) / 2, NA_real_)

group_disp <- data.frame(
  school       = levels(groups),
  n_j          = n_j,
  b_intercept  = round(b_j[, re_names[1L]], 3),
  sigma2_hat_j = round(sigma2_hat_j, 3),
  ratio_pooled = round(sigma2_hat_j / sigma2_pooled, 3),
  win_shape_j  = round(prior_shape_j, 3),
  win_rate_j   = round(prior_rate_j, 1),
  row.names    = NULL
)

cat(sprintf(
  "\n=== Exploratory: implied per-school dispersion (all %d full-rank schools) ===\n\n",
  nlevels(groups)
))
cat(sprintf("  Pooled dispersion_ranef (sigma(lmer)^2): %.4f\n\n", sigma2_pooled))
print(group_disp)
cat(sprintf(
  paste0(
    "\n  win_shape_j/win_rate_j: per-school ING window Gamma(shape, rate) on\n",
    "  1/sigma2_j, mean-matched to THAT school's own sigma2_hat_j with shape\n",
    "  set by its own df_j = n_j - %d (win_rate_j/(win_shape_j-1) == sigma2_hat_j).\n",
    "  This is a per-school reference/window quantity, not the actual dGamma()\n",
    "  prior fed to the sampler -- see the SHARED (a0, b0) hyperprior below.\n"
  ),
  p_re
))

valid <- !is.na(sigma2_hat_j)
cat(sprintf(
  "\n  Per-school sigma2_hat_j (n = %d of %d schools with n_j > p_re = %d):\n",
  sum(valid), nlevels(groups), p_re
))
cat(sprintf(
  "    min = %.3f   median = %.3f   mean = %.3f   max = %.3f   sd = %.3f\n",
  min(sigma2_hat_j[valid]), stats::median(sigma2_hat_j[valid]),
  mean(sigma2_hat_j[valid]), max(sigma2_hat_j[valid]), stats::sd(sigma2_hat_j[valid])
))
cat(sprintf(
  "  pooled = %.3f;  per-school ratio range = [%.3f, %.3f]\n\n",
  sigma2_pooled,
  min(sigma2_hat_j[valid]) / sigma2_pooled,
  max(sigma2_hat_j[valid]) / sigma2_pooled
))

## ---------------------------------------------------------------------------
## EXPLORATORY: empirical-Bayes hyperprior for sigma2_j from the per-school
## classical estimates above.
##
## Model:  1/sigma2_j ~ Gamma(a0, rate = b0) iid across schools j = 1..k
##         (same Inverse-Gamma family dGamma() already uses)
##         RSS_j | sigma2_j ~ sigma2_j * ChiSq(df_j),  df_j = n_j - p_re
##
## mu_hat = E[sigma2_j] is exactly the pooled dispersion_ranef already
## computed above (the RSS-weighted mean -- see derivation below).
## tau2_hat = Var[sigma2_j] (the between-school heterogeneity we actually
## want) is estimated with a DerSimonian-Laird-style moment correction: it
## subtracts the *expected within-school sampling variance*
## 2*sigma2_j^2/df_j (exact for a chi-square estimator) from the observed
## spread of sigma2_hat_j, attributing only the remainder to real
## between-school heterogeneity.  This is the same estimator used for tau^2
## in random-effects meta-analysis, applied here to variances instead of
## means.
##
## Weights w_j = df_j / (2*mu_hat^2) are proportional to 1/Var(sigma2_hat_j);
## note sum(w_j * sigma2_hat_j) / sum(w_j) collapses exactly to
## sum(RSS_j)/sum(df_j) = mu_hat, i.e. the already-computed pooled estimate
## IS the correct precision-weighted mean for this model -- no separate mean
## estimation step is needed, only the second moment (tau2_hat) is new.
##
## (a0, b0) are then the mean-matched Inverse-Gamma hyperparameters implied
## by (mu_hat, tau2_hat) -- same mean-matching convention already used by
## ing_prior / ing_prior_measurement in Prior_Setup_lmebayes.
## ---------------------------------------------------------------------------
df_j      <- (n_j - p_re)[valid]
sig_valid <- sigma2_hat_j[valid]
k_valid   <- sum(valid)
mu_hat    <- sigma2_pooled

w_j  <- df_j / (2 * mu_hat^2)
Q    <- sum(w_j * (sig_valid - mu_hat)^2)
c_dl <- sum(w_j) - sum(w_j^2) / sum(w_j)
tau2_hat <- max(0, (Q - (k_valid - 1)) / c_dl)

cat(sprintf("\n=== Exploratory: empirical-Bayes hyperprior for sigma2_j ===\n\n"))
cat(sprintf(
  "  DerSimonian-Laird-style heterogeneity test:\n    Q = %.2f  (df = %d under homogeneity; ",
  Q, k_valid - 1
))
cat(sprintf("Q >> df indicates real between-school heterogeneity)\n"))
cat(sprintf("    tau2_hat (between-school Var[sigma2_j]) = %.3f\n", tau2_hat))

if (tau2_hat > 0) {
  a0 <- 2 + mu_hat^2 / tau2_hat
  b0 <- mu_hat * (a0 - 1)
  n0_hyper <- 2 * a0 - 1 - p_re    # same units as n_prior_measurement

  cat(sprintf(
    "\n  Empirical-Bayes hyperprior: 1/sigma2_j ~ Gamma(shape = %.4g, rate = %.4g)\n",
    a0, b0
  ))
  cat(sprintf(
    "    Inverse-Gamma on sigma2_j: mean = %.3f [matches mu_hat], sd = %.3f\n",
    b0 / (a0 - 1), sqrt(b0^2 / ((a0 - 1)^2 * (a0 - 2)))
  ))
  cat(sprintf(
    "    Equivalent n_prior_hyper (same units as n_prior_measurement) = %.3g\n\n",
    n0_hyper
  ))
} else {
  cat("\n  tau2_hat <= 0: no detectable between-school heterogeneity beyond sampling noise.\n\n")
}
## ---------------------------------------------------------------------------

N_DRAWS      <- 50L
K_SCHOOLS    <- c(10L, 15L, 20L, 25L)
PWT_GRID     <- c(0.49, 0.20, 0.01)

results <- list()

for (k in K_SCHOOLS) {
  cat(sprintf(
    "\n\n====================================================\n"
  ))
  cat(sprintf("  Fixture: %d schools\n", k))
  cat(sprintf(
    "====================================================\n\n"
  ))

  fx <- tryCatch(
    .prepare_small5_lmerb_manual(n_schools = k),
    error = function(e) {
      cat("  SKIP fixture k =", k, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fx)) next

  for (pwt in PWT_GRID) {
    cat(sprintf("--- k = %d  pwt_measurement = %.2f ---\n", k, pwt))

    result <- tryCatch({
      ps <- Prior_Setup_lmebayes(
        fx$form,
        data            = fx$dat,
        pwt             = 0.01,
        pwt_measurement = pwt
      )
      pf <- pfamily_list(ps)
      m  <- ps$ing_prior_measurement

      cat(sprintf(
        "  window [%.1f, %.1f]  ratio=%.2f  n_prior=%.4g  n_combined=%.4g\n",
        m$disp_lower, m$disp_upper,
        m$disp_upper / m$disp_lower,
        m$n_prior, m$n_combined
      ))

      disp_pf <- dGamma(
        shape          = m$shape,
        rate           = m$rate,
        beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
        Inv_Dispersion = TRUE,
        disp_lower     = m$disp_lower,
        disp_upper     = m$disp_upper
      )

      elapsed <- system.time({
        fit <- lmerb(
          fx$form,
          data             = fx$dat,
          pfamily_list     = pf,
          dispersion_ranef = disp_pf,
          n                = N_DRAWS,
          progbar          = TRUE,
          verbose          = FALSE
        )
      })["elapsed"]

      cat(sprintf(
        "  elapsed = %.1f s  ranef.iters.mean = %.1f  sigma2.mean = %.1f\n",
        elapsed,
        fit$ranef.iters.mean,
        fit$sigma2.mean
      ))

      list(k = k, pwt = pwt, elapsed = elapsed,
           ranef.iters.mean = fit$ranef.iters.mean,
           sigma2.mean      = fit$sigma2.mean,
           window_lo        = m$disp_lower,
           window_hi        = m$disp_upper,
           ok               = TRUE)

    }, error = function(e) {
      cat(sprintf("  ERROR k=%d pwt=%.2f: %s\n", k, pwt, conditionMessage(e)))
      list(k = k, pwt = pwt, ok = FALSE, error = conditionMessage(e))
    })

    key <- paste0("k", k, "_pwt", gsub("\\.", "", format(pwt, nsmall = 2)))
    results[[key]] <- result
  }
}

## Summary table
cat("\n\n========== SUMMARY ==========\n\n")
cat(sprintf("%-6s %-5s %-10s %-9s %-9s %-16s %s\n",
            "k", "pwt", "elapsed(s)", "cand/draw", "sigma2", "window", "ok"))
for (r in results) {
  if (isTRUE(r$ok)) {
    cat(sprintf("%-6d %-5.2f %-10.1f %-9.1f %-9.1f [%5.1f,%5.1f]  OK\n",
                r$k, r$pwt, r$elapsed, r$ranef.iters.mean, r$sigma2.mean,
                r$window_lo, r$window_hi))
  } else {
    cat(sprintf("%-6d %-5.2f %-10s %-9s %-9s %-16s ERROR: %s\n",
                r$k, r$pwt, "-", "-", "-", "-", r$error))
  }
}

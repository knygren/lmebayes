## Demo: lmerb() full workflow on lme4::sleepstudy (Gaussian)
##
## Small model parallel to example(lmer) in lme4: reaction time versus days
## of sleep deprivation with subject-specific intercepts and slopes.  The
## classical ?lmer example uses correlated RE (Days | Subject); lmebayes
## requires uncorrelated terms, so Days is centered and the formula uses
## (1 + Days_c || Subject), as in data-raw/test_lmer_sleepstudy.R.
##
## Replica of a full lmerb run with stored draws (the ?lmerb man-page example
## uses simulate = FALSE for R CMD check).  Keeps model_setup(),
## Prior_Setup_lmebayes(), lmerb(), print_coef_means(), and short Gaussian
## draws-vs-ICM z diagnostics.
##
##   demo("Ex_14_lmerb_Sleepstudy", package = "lmebayes")

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This demo requires the 'lme4' package.", call. = FALSE)
}

## lmerb: subject random effects on sleepstudy
##
## Gaussian LMM for Reaction (ms).  Fixed effect Days_c (centered); random
## intercept and Days_c slope by Subject.  No level-2 covariates.
##
## Workflow: model_setup(), Prior_Setup_lmebayes(), pfamily_list(), then
## lmerb(pfamily_list = , dispersion_ranef = ).

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)

form <- Reaction ~ Days_c + (1 + Days_c || Subject)

design <- model_setup(form, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
cat("\n=== Prior_Setup_lmebayes ===\n\n")
print(ps)

fit <- lmerb(
  form,
  data             = dat,
  pfamily_list     = pfamily_list(ps),
  dispersion_ranef = ps$dispersion_ranef,
  n                = 1000L
)

cat("\n=== summary(fit) ===\n\n")
print(summary(fit))

lmebayes:::print_coef_means(fit)
cat("\nlme4::fixef(fit$lmer):\n")
print(lme4::fixef(fit$lmer))

re_names <- fit$model_setup$re_coef_names
grp_col  <- fit$model_setup$group_name
grp_levs <- rownames(coef(fit$lmer)[[grp_col]])
n_draws  <- nrow(fit$fixef[[re_names[1L]]])

## --- Block 2 fixed effects: MCMC mean vs ICM posterior mean ------------------
cat("\n=== Posterior means of fixed effects (from Block 2 draws) ===\n\n")
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            "RE component", "parameter", "draws mean", "draws SD", "SE(mean)",
            "ICM mean", "z"))
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            strrep("-", 18L), strrep("-", 28L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L),
            strrep("-", 10L), strrep("-", 7L)))
for (k in re_names) {
  dm_k  <- fit$fixef.means[[k]]
  sd_k  <- apply(fit$fixef[[k]], 2L, sd)
  se_k  <- sd_k / sqrt(n_draws)
  icm_k <- fit$fixef.mode[[k]]
  for (nm in names(dm_k)) {
    z_val <- (dm_k[[nm]] - icm_k[[nm]]) / se_k[[nm]]
    flag  <- if (abs(z_val) > 2) " *" else "  "
    cat(sprintf("  %-18s  %-28s  %10.4f  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                k, nm, dm_k[[nm]], sd_k[[nm]], se_k[[nm]], icm_k[[nm]], z_val, flag))
  }
}
cat("  (* |z| > 2: draws mean inconsistent with exact ICM posterior mean)\n\n")

## --- Random effects: MCMC mean vs ICM posterior mean -----------------------
cat("=== Random effects: MCMC mean vs ICM posterior mean ===\n\n")
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

cat(sprintf("  %-8s  %-14s  %10s  %10s  %10s  %6s\n",
            "group", "RE component", "MCMC mean", "ICM mean", "SE(mean)", "z"))
cat(sprintf("  %-8s  %-14s  %10s  %10s  %10s  %6s\n",
            strrep("-", 8L), strrep("-", 14L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L), strrep("-", 6L)))

n_flagged <- 0L
for (lev in grp_levs) {
  lev_chr <- as.character(lev)
  for (k in re_names) {
    mcmc_m <- re_draws_mean[[lev_chr]][[k]]
    mcmc_s <- re_draws_sd[[lev_chr]][[k]]
    icm_m  <- fit$ranef.mode[lev_chr, k]
    se_val <- mcmc_s / sqrt(n_draws)
    z_val  <- (mcmc_m - icm_m) / se_val
    flag   <- if (abs(z_val) > 3) " *" else "  "
    if (abs(z_val) > 3) n_flagged <- n_flagged + 1L
    cat(sprintf("  %-8s  %-14s  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                lev_chr, k, mcmc_m, icm_m, se_val, z_val, flag))
  }
}
total_tests <- length(grp_levs) * length(re_names)
cat(sprintf(
  "\n  %d of %d tests flagged |z| > 3  (expected ~%.1f by chance at 0.3%% level)\n",
  n_flagged, total_tests, total_tests * 0.003
))
cat("  (* |z| > 3: MCMC mean inconsistent with exact ICM posterior mean)\n")

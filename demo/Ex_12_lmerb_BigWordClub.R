## Demo: lmerb() full workflow on bayesrules::big_word_club
##
## Full model with cross-level moderation (not the simplified ?lmerb man-page
## example).  Preserved as a demo because the run (1000 draws plus factor-level
## diagnostics) takes on the order of a minute.  This demo keeps the complete
## workflow: model_setup(), Prior_Setup_GLMM(), pfamily_list(), lmerb(),
## then draws-vs-ICM z-tests and lmer/mu_all/lmerb factor-level comparisons.
##
## For the lme4-style small model with simulation, see
## demo("Ex_14_lmerb_Sleepstudy", package = "lmebayes").
##
##   demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")

for (pkg in c("bayesrules", "coda")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("This demo requires the '%s' package.", pkg), call. = FALSE)
  }
}

## lmerb: school random effects on bayesrules::big_word_club
##
## Compare factor-level means: lmer coef(), Block 1 prior mean mu_all,
## and posterior means of lmerb draws (averaged over draw within each level).
##
## Student-level random slopes: distracted_a1 and distracted_ppvt (distraction
## during assessment 1 vs during PPVT). Both have fixed main effects + school
## random slopes. Cross-level moderation: free_reduced_lunch:distracted_a1
## (school lunch status moderates the distracted_a1 random slope).
##
## Workflow: model_setup(), Prior_Setup_GLMM(), pfamily_list(), then
## lmerb(pfamily_list = , dispersion_ranef = ).

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)


form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)



design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_GLMM(form_lmer, data = dat, pwt = 0.01)
cat("\n=== Prior_Setup_GLMM ===\n\n")
print(ps)

## Defaults: tv_tol = 0.01 (each stored draw within 0.01 TV of the exact
## joint posterior; m_convergence derived from the Nygren (2020) Theorem 3
## bound), no fixed seed.
fit <- lmerb(
  form_lmer,
  data = dat,
  pfamily_list = pfamily_list(ps),
  dispersion_ranef = ps$dispersion_ranef,
  n = 1000L
)

cat("\n=== summary(fit) ===\n\n")
print(summary(fit))

grp_col  <- fit$model_setup$group_name
re_names <- fit$model_setup$re_coef_names
grp_levs <- rownames(coef(fit$lmer)[[grp_col]])

## --- Posterior means of fixed effects from Block 2 draws -----------------
## fit$fixef.means[[k]]: MCMC mean of the n Block 2 gamma_k draws.
## Directly comparable to lme4::fixef() (same parameter scale).
## fit$fixef.mode: exact ICM posterior mean (used as H0 for z-test below).
## z = (draws mean - ICM mean) / (draws SD / sqrt(n)); flag |z| > 2 with *.
cat("\n=== Posterior means of fixed effects (from Block 2 draws) ===\n\n")
n_draws     <- nrow(fit$fixef[[re_names[1L]]])
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            "RE component", "parameter", "draws mean", "draws SD", "SE(mean)", "ICM mean", "z"))
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
cat("  (* |z| > 2: draws mean inconsistent with exact ICM posterior mean)\n")
cat("\n")

## --- Z-test: MCMC mean vs ICM posterior mean for random effects -------------
## fit$ranef.mode: exact ICM posterior mean (J x p_re), rows = group levels.
## MCMC mean computed as the per-group average of the n Block 1 draws.
## z = (MCMC mean - ICM mean) / (draws SD / sqrt(n))
## With J * p_re tests we expect a few |z| > 2 by chance; flag |z| > 3.

cat("=== Random effects: MCMC mean vs ICM posterior mean ===\n\n")

## Compute per-group MCMC means and SDs from the stored coefficients.
## simplify=FALSE forces a named list regardless of output length.
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
n_re_draws <- n_draws   # same n as fixed effects

cat(sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
            "group", "RE component", "MCMC mean", "ICM mean", "SE(mean)", "z"))
cat(sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
            strrep("-", 6L), strrep("-", 18L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L), strrep("-", 6L)))

icm_b <- fit$ranef.mode   # J x p_re, rownames = group levels

n_flagged <- 0L
for (lev in grp_levs) {
  lev_chr <- as.character(lev)
  for (k in re_names) {
    mcmc_m <- re_draws_mean[[lev_chr]][[k]]
    mcmc_s <- re_draws_sd[[lev_chr]][[k]]
    icm_m  <- icm_b[lev_chr, k]
    se_val <- mcmc_s / sqrt(n_re_draws)
    z_val  <- (mcmc_m - icm_m) / se_val
    flag   <- if (abs(z_val) > 3) " *" else "  "
    if (abs(z_val) > 3) n_flagged <- n_flagged + 1L
    cat(sprintf("  %-6s  %-18s  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                lev_chr, k, mcmc_m, icm_m, se_val, z_val, flag))
  }
}
total_tests <- length(grp_levs) * length(re_names)
cat(sprintf(
  "\n  %d of %d tests flagged |z| > 3  (expected ~%.1f by chance at 0.3%% level)\n",
  n_flagged, total_tests, total_tests * 0.003
))
cat("  (* |z| > 3: MCMC mean inconsistent with exact ICM posterior mean)\n\n")

cat("\n=== Random effects: lmer reference vs lmerb chain mean ===\n\n")
lmebayes:::print_mer_bayes_re_compare(fit, detail = "full", digits = 3L)

## coda: one chain per factor level (10000 draws x p_re per school)
coef_split <- split(
  fit$coefficients,
  factor(fit$coefficients[[grp_col]], levels = grp_levs)
)
mcmc_by_level <- coda::mcmc.list(lapply(coef_split, function(d) {
  coda::mcmc(as.matrix(d[, re_names, drop = FALSE]))
}))
ex_level <- as.character(grp_levs[1L])
cat("\ncoda summary for factor level", ex_level, ":\n")
print(summary(mcmc_by_level[[ex_level]]))

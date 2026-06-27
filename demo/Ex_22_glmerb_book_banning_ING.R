## Demo: glmerb() binomial workflow on bayesrules::book_banning (ING priors)
##
## Same model as demo/Ex_16_glmerb_book_banning.R, but Block~2 RE components use
## dIndependent_Normal_Gamma (sampled tau^2_k) instead of dNormal (fixed dispersion).
## Non-Gaussian ING uses the R sweep-outer driver (run_sweep_outer_chains_v6):
## pilot stage, then main chains from pilot colMeans + tau^2 plug-in.
##
## After fitting:
##   print(fit, sweep_history = TRUE, max_sweeps = 5)
##   plot_sweep_history_diag(fit$sweep_history$main, coef_focus, what = "mean")
##
##   demo("Ex_22_glmerb_book_banning_ING", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

data(book_banning, package = "bayesrules", envir = environment())

dat <- book_banning[, c(
  "state", "removed", "violent", "antifamily", "language"
)]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)

form_book <- removed ~ violent + antifamily + language + (1 | state)
form_glmerb <- removed_i ~ violent_i + (1 + violent_i || state)

ps <- Prior_Setup_lmebayes(
  form_glmerb,
  data           = dat,
  family         = binomial(),
  pwt            = 0.01,
  pwt_dispersion = 0.2
)
cat("\n=== Prior_Setup_lmebayes (ING calibration) ===\n\n")
print(ps)

pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

set.seed(42L)
fit <- glmerb(
  form_glmerb,
  data         = dat,
  family       = binomial(),
  pfamily_list = pf,
  n            = 3000L,
  mode_gap_max = 1.0,
  progbar      = TRUE
)

stopifnot(isTRUE(fit$prior$any_non_normal))
stopifnot(is.matrix(fit$fixef.dispersion))
stopifnot(all(is.finite(fit$fixef.dispersion)), all(fit$fixef.dispersion > 0))

cat("m_convergence used:", fit$convergence$m_convergence, "\n")
cat(sprintf(
  "Pilot vs mode (chi-squared): p = %.4g (n_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot
))

cat("\n--- Sampled RE variance (tau^2) posterior means ---\n")
print(fit$fixef.dispersion.mean)

cat("\n--- Ch. 18 reference glmer (random intercept; all three reasons) ---\n")
fit_book <- lme4::glmer(form_book, data = dat, family = binomial())
print(lme4::fixef(fit_book))

lmebayes:::print_coef_means(fit)
print(fit)
summary(fit)

coef_focus <- list(
  c("(Intercept)", "(Intercept)"),
  c("violent_i", "(Intercept)")
)

for (st in list(fit$sweep_history$pilot, fit$sweep_history$main)) {
  if (is.null(st)) next
  plot_sweep_history_diag(st, coef_focus)
}

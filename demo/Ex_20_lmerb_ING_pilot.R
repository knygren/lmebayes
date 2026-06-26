## Demo: lmerb() with dIndependent_Normal_Gamma Block~2 priors + pilot stage
##
## Gaussian LMM with ING random-effect dispersion: runs the two-stage workflow
## (pilot chains from ICM start, main chains from pilot mean) like glmerb().
## Use this as a regression target when changing ING start points or sampling.
##
##   demo("Ex_20_lmerb_ING_pilot", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This demo requires the 'lme4' package.", call. = FALSE)
}

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(dat, !is.na(score_ppvt))

## School random intercept; one ING component (same structure as test_ing_sampling.R).
form <- score_ppvt ~ private_school + (1 | school_id)

ps <- Prior_Setup_lmebayes(
  form,
  data             = dat,
  pwt              = 0.01,
  pwt_dispersion   = 0.2
)
cat("\n=== Prior_Setup_lmebayes (ING calibration) ===\n\n")
print(ps)

pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

cat("\n=== lmer reference fit ===\n\n")
fit_lmer <- lme4::lmer(form, data = dat, REML = TRUE)
print(summary(fit_lmer))

## gap_tol = 0.05 => n_pilot ~ 16 (legacy Hotelling bound); main n = 500.
fit <- lmerb(
  form,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = ps$dispersion_ranef,
  n                = 10000L,
  gap_tol          = 0.05,
  mode_gap_max     = 1.0
)

stopifnot(isTRUE(fit$prior$any_non_normal))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(fit$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))
stopifnot(is.finite(fit$pilot_chisq$p_value))

cat("\n=== summary(lmerb fit) ===\n\n")
print(summary(fit))

stopifnot(!is.null(fit$sweep_history))
stopifnot(!is.null(fit$sweep_history$pilot))
stopifnot(!is.null(fit$sweep_history$main))

cat("\n=== Block~2 sweep summaries (pilot, then main) ===\n\n")
print(fit$sweep_history$pilot)
print(fit$sweep_history$main)

re_names <- fit$model_setup$re_coef_names
pr_int   <- pf[["(Intercept)"]]$prior_list
t2_draws <- fit$fixef.dispersion[, "(Intercept)"]
stopifnot(
  all(is.finite(t2_draws)), all(t2_draws > 0),
  all(t2_draws >= pr_int$disp_lower),
  all(t2_draws <= pr_int$disp_upper),
  stats::sd(t2_draws) > 0
)

cat(sprintf(
  "\nPilot vs mode (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  fit$convergence$m_convergence_pilot
))
cat(sprintf(
  "tau^2 posterior mean = %.4f  [truncation window (%.4f, %.4f)]\n",
  fit$fixef.dispersion.mean[["(Intercept)"]],
  pr_int$disp_lower,
  pr_int$disp_upper
))

## Centering: main-stage mean vs pilot mean vs plug-in start (Block~2 hyperparameters).
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
beta_bar <- unlist(lapply(re_names, function(k) fit$fixef.means[[k]]))
names(beta_bar) <- cn
theta_mode  <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
theta_pilot <- unlist(lapply(re_names, function(k) fit$fixef.init[[k]]))
names(theta_mode) <- names(theta_pilot) <- cn

X <- do.call(cbind, lapply(re_names, function(k) fit$fixef[[k]]))
colnames(X) <- cn
V <- stats::cov(X)
V_inv <- tryCatch(
  solve(V),
  error = function(e) solve(V + diag(1e-8 * mean(diag(V)), ncol(V)))
)
d_pilot <- beta_bar - theta_pilot
d_mode  <- beta_bar - theta_mode
Q_pilot <- as.numeric(t(d_pilot) %*% V_inv %*% d_pilot)
Q_mode  <- as.numeric(t(d_mode) %*% V_inv %*% d_mode)
p_pilot <- stats::pchisq(Q_pilot, df = ncol(X), lower.tail = FALSE)
p_mode  <- stats::pchisq(Q_mode, df = ncol(X), lower.tail = FALSE)

cat(sprintf(
  "Overall centering (chi-squared): p(mean=pilot)=%.4g, p(mean=mode)=%.4g\n",
  p_pilot, p_mode
))

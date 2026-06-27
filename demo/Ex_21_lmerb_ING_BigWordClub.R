## Demo: lmerb() with ING Block~2 priors + pilot on the full big_word_club model
##
## Combines demo/Ex_12_lmerb_BigWordClub.R (school random intercept and slopes,
## cross-level moderation) with demo/Ex_20_lmerb_ING_pilot.R (dIndependent_Normal_
## Gamma dispersion, two-stage pilot/main sampling).
##
## ING sampling (non-dNormal Block~2) uses glmbayesCore's R sweep-outer driver
## run_sweep_outer_chains_v6 (pilot, then main). diag_sweeps = TRUE auto-prints
## one combined Block~2 chain-mean table per stage; sweep_history is always
## stored (print(fit$sweep_history$pilot) / fit$sweep_history$main).
##
##   demo("Ex_21_lmerb_ING_BigWordClub", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This demo requires the 'lme4' package.", call. = FALSE)
}

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

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data             = dat,
  pwt              = 0.01,
  pwt_dispersion   = 0.2
)
cat("\n=== Prior_Setup_lmebayes (ING calibration) ===\n\n")
print(ps)

pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

cat("\n=== lmer reference fit ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
print(summary(fit_lmer))

## diag_sweeps = TRUE: auto-print per-stage sweep tables (progbar defaults FALSE).
## gap_tol = 0.05 => n_pilot from the Hotelling bound (~16 for this model).
fit <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = ps$dispersion_ranef,
  n                = 3000L,
  gap_tol          = 0.05,
  mode_gap_max     = 1.0,
  diag_sweeps      = FALSE
)

stopifnot(isTRUE(fit$prior$any_non_normal))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(!is.null(fit$sweep_history))
stopifnot(!is.null(fit$sweep_history$pilot))
stopifnot(!is.null(fit$sweep_history$main))
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

## tau^2_k per RE component: sampled inside [disp_lower, disp_upper].
for (k in re_names) {
  pr_k <- pf[[k]]$prior_list
  t2   <- fit$fixef.dispersion[, k]
  stopifnot(
    all(is.finite(t2)), all(t2 > 0),
    all(t2 >= pr_k$disp_lower),
    all(t2 <= pr_k$disp_upper),
    stats::sd(t2) > 0
  )
  cat(sprintf(
    "\n%s tau^2: post mean = %.4f  [window (%.4f, %.4f); plugin disp_lower = %.4f]\n",
    k,
    fit$fixef.dispersion.mean[[k]],
    pr_k$disp_lower,
    pr_k$disp_upper,
    pr_k$disp_lower
  ))
}

cat(sprintf(
  "\nPilot vs plug-in start (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  fit$convergence$m_convergence_pilot
))

## Block~2 hyperparameters: prior mean, gamma @ lmer tau2, pilot init, MCMC mean.
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
beta_bar <- unlist(lapply(re_names, function(k) fit$fixef.means[[k]]))
theta_plug <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
theta_prior <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  unname(fit$prior$prior_list[[k]]$mu_fixef[nms])
}))
theta_pilot <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  unname(fit$fixef.init[[k]][nms])
}))
names(beta_bar) <- names(theta_plug) <- names(theta_prior) <- names(theta_pilot) <- cn

block2_cmp <- data.frame(
  prior_mean      = unname(theta_prior),
  gamma_lmer_tau2 = unname(theta_plug),
  pilot_mean      = unname(theta_pilot),
  mcmc_mean       = unname(beta_bar),
  row.names       = cn,
  check.names     = FALSE
)
cat("\n=== Block 2 hyperparameters (prior / plug-in / pilot / MCMC) ===\n\n")
print(round(block2_cmp, 4))

## Multivariate centering: main mean vs pilot start vs plug-in start.
X <- do.call(cbind, lapply(re_names, function(k) fit$fixef[[k]]))
colnames(X) <- cn
V <- stats::cov(X)
V_inv <- tryCatch(
  solve(V),
  error = function(e) solve(V + diag(1e-8 * mean(diag(V)), ncol(V)))
)
d_pilot <- beta_bar - theta_pilot
d_plug  <- beta_bar - theta_plug
Q_pilot <- as.numeric(t(d_pilot) %*% V_inv %*% d_pilot)
Q_plug  <- as.numeric(t(d_plug) %*% V_inv %*% d_plug)
p_pilot <- stats::pchisq(Q_pilot, df = ncol(X), lower.tail = FALSE)
p_plug  <- stats::pchisq(Q_plug, df = ncol(X), lower.tail = FALSE)

cat(sprintf(
  "\nOverall centering (chi-squared, p = %d hyperparameters): p(mean=pilot)=%.4g, p(mean=plug-in)=%.4g\n",
  ncol(X), p_pilot, p_plug
))

pilot_mean <- unname(theta_pilot)
post_mean  <- unname(beta_bar)
post_sd    <- unlist(lapply(re_names, function(k) apply(fit$fixef[[k]], 2L, sd)))
plug_in    <- unname(theta_plug)
names(post_sd) <- cn

n_main <- nrow(fit$fixef[[re_names[1L]]])
mc_se  <- post_sd / sqrt(n_main)

tab <- data.frame(
  pilot_mean  = round(pilot_mean, 4),
  post_mean   = round(post_mean, 4),
  difference  = round(post_mean - pilot_mean, 4),
  post_sd     = round(post_sd, 4),
  mc_se       = round(mc_se, 4),
  z_vs_pilot  = round((post_mean - pilot_mean) / mc_se, 2),
  plug_in     = round(plug_in, 4),
  row.names   = cn,
  check.names = FALSE
)

cat("\n=== Block 2: pilot mean vs posterior mean ===\n\n")
print(tab)

## Sweep-history diagnostics: cross-chain mean and SD vs inner sweep (pilot and main).
## Helps spot coefficients whose chains spread or drift across inner sweeps.
coef_focus <- list(
  c("(Intercept)", "(Intercept)"),
  c("(Intercept)", "private_school"),
  c("(Intercept)", "title1"),
  c("(Intercept)", "free_reduced_lunch"),
  c("distracted_ppvt", "(Intercept)"),
  c("distracted_a1", "(Intercept)"),
  c("distracted_a1", "free_reduced_lunch")
)

for (st in list(fit$sweep_history$pilot, fit$sweep_history$main)) {
  if (is.null(st)) next
  plot_sweep_history_diag(st, coef_focus)
}

# Manual validation: lmerb dGamma sigma^2 — Block~2 fixef + RE (mer_full vs chain mean).
# Run after Block~1 / BlockEnvelopeCentering / dGamma sampling changes.
#
# Same Big Word Club model as test_lmerb_mer_re_validation.R; only the
# measurement-dispersion route differs (dGamma sigma^2 vs fixed scalar).
#
# Sampler: lmerb() defaults (tv_tol, gap_tol, mode_gap_max). No set.seed().
# Main draws: n = 1000 (dGamma, fixed tau^2); n = 3000 (dGamma + ING tau^2).
# progbar = TRUE only.
#
# Routes:
#   rLMMindepNormalGamma_reg_known_vcov      (parallel to Gaussian §1)
#   rLMMindepNormalGamma_reg_estimated_vcov   (parallel to ING §2)
#
#   Rscript tests/manual/test_lmerb_dgamma_mer_re_validation.R

source("tests/manual/_load.R")
source("tests/manual/_lmerb_dgamma_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_DGAMMA      <- 1000L
N_DGAMMA_ING  <- 3000L

## --- Shared Big Word Club model (see _bwc_lmerb_fixture.R) -------------------

fx <- .prepare_bwc_lmerb_manual(pwt = 0.01)
dat  <- fx$dat
form <- fx$form
expected_re <- fx$design$re_coef_names
stopifnot(identical(
  expected_re,
  c("(Intercept)", "distracted_ppvt", "distracted_a1")
))

## --- 1. Random sigma^2, fixed Block~2 tau^2 (parallel to Gaussian §1) -------

pf <- pfamily_list(fx$ps)
disp_pf <- .dgamma_dispersion_ranef(fx$ps, narrow_window = TRUE)

cat("\n=== lmerb dGamma sigma^2, fixed Block~2 tau^2; n =", N_DGAMMA, "===\n\n")
fit <- lmerb(
  form,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = disp_pf,
  n                = N_DGAMMA,
  progbar          = TRUE
)

stopifnot(inherits(fit, "lmerb"))
stopifnot(identical(fit$prior$dispersion_mode, "gamma"))
stopifnot(!isTRUE(fit$prior$any_non_normal))
stopifnot(is.null(fit$fixef.dispersion))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(fit$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))
stopifnot(is.finite(fit$pilot_chisq$p_value))
stopifnot(!is.null(fit$sweep_history$main))
re_names <- fit$model_setup$re_coef_names
stopifnot(identical(re_names, expected_re))
stopifnot(identical(nrow(fit$fixef[[re_names[1L]]]), N_DGAMMA))

cat(sprintf(
  "\nPilot vs mode: p = %.4g (n_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot
))

lmebayes:::.validate_lmerb_block2_fixef_lmer(fit, label = "dGamma fixed tau^2")

.validate_lmerb_re(fit, label = "dGamma fixed tau^2", cor_full_min = 0.85)

## --- 2. Random sigma^2 + ING Block~2 (parallel to ING §2) --------------------

ps2 <- Prior_Setup_lmebayes(
  form,
  data             = dat,
  pwt              = 0.01,
  pwt_dispersion   = 0.2
)
pf2 <- pfamily_list(ps2, ptypes = "dIndependent_Normal_Gamma")
disp_pf2 <- .dgamma_dispersion_ranef(ps2, narrow_window = FALSE)

cat("\n=== lmerb dGamma sigma^2 + ING Block~2; n =", N_DGAMMA_ING, "===\n\n")
fit2 <- lmerb(
  form,
  data             = dat,
  pfamily_list     = pf2,
  dispersion_ranef = disp_pf2,
  n                = N_DGAMMA_ING,
  progbar          = TRUE
)

stopifnot(identical(fit2$prior$dispersion_mode, "gamma"))
stopifnot(isTRUE(fit2$prior$any_non_normal))
stopifnot(!is.null(fit2$fixef.dispersion))
stopifnot(!is.null(fit2$pilot_chisq))
stopifnot(fit2$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit2$pilot_chisq$n_pilot, fit2$convergence$n_pilot))
stopifnot(is.finite(fit2$pilot_chisq$p_value))
stopifnot(!is.null(fit2$sweep_history$pilot))
stopifnot(!is.null(fit2$sweep_history$main))
re_names2 <- fit2$model_setup$re_coef_names
stopifnot(identical(re_names2, expected_re))
stopifnot(identical(re_names2, re_names))
stopifnot(identical(nrow(fit2$fixef[[re_names2[1L]]]), N_DGAMMA_ING))

for (k in re_names2) {
  pr_k <- pf2[[k]]$prior_list
  t2   <- fit2$fixef.dispersion[, k]
  stopifnot(
    all(is.finite(t2)), all(t2 > 0),
    all(t2 >= pr_k$disp_lower),
    all(t2 <= pr_k$disp_upper),
    stats::sd(t2) > 0
  )
}

cat(sprintf(
  "\nPilot vs plug-in start: p = %.4g (n_pilot = %d)\n",
  fit2$pilot_chisq$p_value,
  fit2$pilot_chisq$n_pilot
))

lmebayes:::.validate_lmerb_block2_fixef_lmer(fit2, label = "dGamma + ING")

.validate_lmerb_re(fit2, label = "dGamma + ING", cor_full_min = 0.85)

cat("\ntest_lmerb_dgamma_mer_re_validation: OK\n")

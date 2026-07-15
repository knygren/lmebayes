# Manual smoke test: dGamma Block~1 on a tiny 5-school model.
#
# Fixed: intercept + distracted_ppvt (Prior_Setup requires fixed main for RE slope).
# RE: (Intercept) + distracted_ppvt (2 columns).  5 schools x ~9 faces (Gridtype 3, l1 = 2).
#
# Requires dev glmbayesCore (sibling package) for pwt_measurement and Block~1 ING.
#   cd .../lmebayes
#   Rscript tests/manual/test_lmerb_dgamma_small5_validation.R

source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_DGAMMA <- 1000L

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
dat  <- fx$dat
form <- fx$form
expected_re <- fx$design$re_coef_names

stopifnot(identical(expected_re, c("(Intercept)", "distracted_ppvt")))
stopifnot(nlevels(fx$design$groups) == 5L)

source("tests/manual/_reference_mer_compare.R")
ref_fits <- .print_reference_mer_compare(
  form        = form,
  dat         = dat,
  group_name  = "school_id",
  dispformula = ~1
)

ps <- Prior_Setup_lmebayes(
  form,
  data            = dat,
  pwt             = 0.01,
  pwt_measurement = 0.49
)
pf <- pfamily_list(ps)

m_disp <- ps$ing_prior_measurement
cat(sprintf(
  "\n=== ING sigma^2 prior (pwt_measurement = %.2f); window [%.1f, %.1f] ===\n\n",
  ps$pwt_measurement,
  m_disp$disp_lower,
  m_disp$disp_upper
))

disp_pf <- dGamma(
  shape          = m_disp$shape,
  rate           = m_disp$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m_disp$disp_lower,
  disp_upper     = m_disp$disp_upper
)

cat("\n=== Small 5-school dGamma smoke test; n =", N_DGAMMA, "===\n\n")
cat("draw_engine after fit should be rGLMM_sweep_ing_block1_ind\n\n")

t_fit <- system.time({
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pf,
    dispersion_ranef = disp_pf,
    n                = N_DGAMMA,
    progbar          = TRUE,
    verbose          = TRUE
  )
})

cat(sprintf("\n=== Timing: lmerb elapsed = %.2f s ===\n", t_fit["elapsed"]))

summary(fit)

stopifnot(inherits(fit, "lmerb"))
stopifnot(identical(fit$prior$dispersion_mode, "gamma"))
stopifnot(identical(fit$draw_engine, "rGLMM_sweep_ing_block1_ind"))

re_names <- fit$model_setup$re_coef_names
stopifnot(identical(re_names, expected_re))
stopifnot(identical(nrow(fit$fixef[[re_names[1L]]]), N_DGAMMA))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(is.finite(fit$pilot_chisq$p_value))

# Block~1 sigma^2 draws (dGamma measurement dispersion)
stopifnot(is.numeric(fit$sigma2), length(fit$sigma2) == N_DGAMMA)
stopifnot(is.finite(fit$sigma2.mean))
cat(sprintf(
  "sigma2: mean = %.4f, sd = %.4f\n",
  fit$sigma2.mean, stats::sd(fit$sigma2)
))

cat(sprintf(
  "\nPilot vs mode: p = %.4g (n_pilot = %d, m_convergence = %d, m_pilot = %s)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  fit$m_convergence,
  if (is.null(fit$m_convergence_pilot)) NA else fit$m_convergence_pilot
))

cat("\ntest_lmerb_dgamma_small5_validation: OK\n")

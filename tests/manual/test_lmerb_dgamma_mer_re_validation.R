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

# Requires dev glmbayesCore (sibling package) for pwt_measurement and Block~1 ING.
#
#   Rscript tests/manual/test_lmerb_dgamma_mer_re_validation.R

source("tests/manual/_load.R")
source("tests/manual/_bwc_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)
source("tests/manual/_block2_fixef_validate.R")
.bind_manual_block2_fixef()
source("tests/manual/_lmerb_re_validate.R")
.bind_lmerb_re_validate()



N_DGAMMA      <- 1000L
N_DGAMMA_ING  <- 3000L
PWT_MEASUREMENT <- 0.2

## Block~1: BlockEnvelopeCentering → Build → DispersionBuild → Sim; ranef.mode / fixef.mode
## are plug-in ICM — skip ICM thresholds §1–§2.
DGAMMA_ENV_RE <- list(cor_full_min = 0.85, cor_icm_min = 0, mode_match_min = 0)



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

ps <- Prior_Setup_GLMM(
  form,
  data            = dat,
  pwt             = 0.01,
  pwt_measurement = PWT_MEASUREMENT
)
pf <- pfamily_list(ps)

m_disp <- ps$ing_prior_measurement
cat(sprintf(
  paste0(
    "\n=== dGamma sigma^2 prior (section 1; pwt_measurement = %.2g); ",
    "window [%.4g, %.4g]; sigma2_hat = %.4g ===\n\n"
  ),
  ps$pwt_measurement,
  m_disp$disp_lower,
  m_disp$disp_upper,
  m_disp$sigma2_hat
))

disp_pf <- dGamma(
  shape          = m_disp$shape,
  rate           = m_disp$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m_disp$disp_lower,
  disp_upper     = m_disp$disp_upper
)

cat("\n=== lmerb dGamma sigma^2, fixed Block~2 tau^2; n =", N_DGAMMA, "===\n\n")



fit <- lmerb(

  form,

  data             = dat,

  pfamily_list     = pf,

  dispersion_ranef = disp_pf,

  n                = N_DGAMMA,

  progbar          = TRUE

)

summary(fit)

stopifnot(inherits(fit, "lmerb"))

stopifnot(identical(fit$prior$dispersion_mode, "gamma"))

stopifnot(!isTRUE(fit$prior$any_non_normal))

re_names <- fit$model_setup$re_coef_names

stopifnot(identical(re_names, expected_re))

stopifnot(is.matrix(fit$fixef.dispersion))

stopifnot(

  identical(nrow(fit$fixef.dispersion), N_DGAMMA),

  identical(colnames(fit$fixef.dispersion), re_names),

  all(is.finite(fit$fixef.dispersion)), all(fit$fixef.dispersion > 0)

)

for (k in re_names) {

  t2 <- fit$fixef.dispersion[, k]

  stopifnot(stats::sd(t2) == 0)

  stopifnot(isTRUE(all.equal(

    unique(t2),

    pf[[k]]$prior_list$dispersion,

    tolerance = 1e-6

  )))

}

stopifnot(!is.null(fit$pilot_chisq))

stopifnot(fit$pilot_chisq$n_pilot > 0L)

stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))

stopifnot(is.finite(fit$pilot_chisq$p_value))

stopifnot(!is.null(fit$sweep_history$main))

stopifnot(identical(nrow(fit$fixef[[re_names[1L]]]), N_DGAMMA))



cat(sprintf(

  "\nPilot vs mode: p = %.4g (n_pilot = %d)\n",

  fit$pilot_chisq$p_value,

  fit$pilot_chisq$n_pilot

))



.validate_manual_block2_fixef(
  fit, label = "dGamma fixed tau^2", z_icm_max = Inf
)

do.call(
  lmebayes:::.validate_lmerb_re,
  c(list(fit = fit, label = "dGamma fixed tau^2"), DGAMMA_ENV_RE)
)



## --- 2. Random sigma^2 + ING Block~2 (parallel to ING §2) --------------------



ps2 <- Prior_Setup_GLMM(

  form,

  data             = dat,

  pwt              = 0.01,

  pwt_dispersion   = 0.2,

  pwt_measurement  = PWT_MEASUREMENT

)

pf2 <- pfamily_list(ps2, ptypes = "dIndependent_Normal_Gamma")

m_disp2 <- ps2$ing_prior_measurement

cat(sprintf(

  paste0(

    "\n=== dGamma sigma^2 prior (section 2; pwt_measurement = %.2g); ",

    "window [%.4g, %.4g]; sigma2_hat = %.4g ===\n\n"

  ),

  ps2$pwt_measurement,

  m_disp2$disp_lower,

  m_disp2$disp_upper,

  m_disp2$sigma2_hat

))

disp_pf2 <- dGamma(

  shape          = m_disp2$shape,

  rate           = m_disp2$rate,

  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),

  Inv_Dispersion = TRUE,

  disp_lower     = m_disp2$disp_lower,

  disp_upper     = m_disp2$disp_upper

)



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



.validate_manual_block2_fixef(
  fit2, label = "dGamma + ING", z_icm_max = Inf
)

do.call(
  lmebayes:::.validate_lmerb_re,
  c(list(fit = fit2, label = "dGamma + ING"), DGAMMA_ENV_RE)
)



cat("\ntest_lmerb_dgamma_mer_re_validation: OK\n")



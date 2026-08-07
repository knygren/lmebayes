# Manual validation: lmerb Gibbs + mer_full vs chain-mean RE.

# Run after major sampler / Block~2 / centering changes.

#

# Sampler: lmerb() defaults (tv_tol, gap_tol, mode_gap_max). No set.seed().

# Main draws: n = 1000 (Gaussian); n = 3000 (ING). progbar = TRUE only.

#

# Routes:

#   rLMMNormal_reg_known_vcov       (Big Word Club; Block~2 vs lmer + RE cor)

#   rLMMNormal_reg_estimated_vcov   (Ex_21-style ING; same model as §1, ING priors)

#

#   Rscript tests/manual/test_lmerb_mer_re_validation.R



source("tests/manual/_load.R")
source("tests/manual/_bwc_lmerb_fixture.R")
.manual_test_load()
source("tests/manual/_block2_fixef_validate.R")
.bind_manual_block2_fixef()
source("tests/manual/_lmerb_re_validate.R")
.bind_lmerb_re_validate()



N_GAUSSIAN <- 1000L

N_ING      <- 3000L



## --- Shared Big Word Club model (see _bwc_lmerb_fixture.R) -------------------



fx <- .prepare_bwc_lmerb_manual(pwt = 0.01)

dat  <- fx$dat

form <- fx$form

expected_re <- fx$design$re_coef_names

stopifnot(identical(

  expected_re,

  c("(Intercept)", "distracted_ppvt", "distracted_a1")

))



## --- 1. Gaussian known vcov (Ex_12-style) ------------------------------------



ps <- fx$ps



cat("\n=== lmerb Gaussian known vcov; n =", N_GAUSSIAN, "===\n\n")

fit <- lmerb(

  form,

  data             = dat,

  pfamily_list     = pfamily_list(ps),

  dispersion_ranef = ps$dispersion_ranef,

  n                = N_GAUSSIAN,

  progbar          = TRUE

)



stopifnot(inherits(fit, "lmerb"))

re_names <- fit$model_setup$re_coef_names

stopifnot(identical(re_names, expected_re))

n_draws <- nrow(fit$fixef[[re_names[1L]]])

stopifnot(identical(n_draws, N_GAUSSIAN))



.validate_manual_block2_fixef(fit, label = "Gaussian known vcov")



.validate_lmerb_re(fit, label = "Gaussian known vcov")



## --- 2. ING estimated vcov (Ex_21-style; same model + data, ING priors) ------



ps2 <- Prior_Setup_GLMM(

  form,

  data             = dat,

  pwt              = 0.01,

  pwt_dispersion   = 0.2

)

pf2 <- pfamily_list(ps2, ptypes = "dIndependent_Normal_Gamma")



cat("\n=== lmerb ING estimated vcov; n =", N_ING, "===\n\n")

fit2 <- lmerb(

  form,

  data             = dat,

  pfamily_list     = pf2,

  dispersion_ranef = ps2$dispersion_ranef,

  n                = N_ING,

  progbar          = TRUE

)



stopifnot(inherits(fit2, "lmerb"))

stopifnot(isTRUE(fit2$prior$any_non_normal))

stopifnot(fit2$pilot_chisq$n_pilot > 0L)

re_names2 <- fit2$model_setup$re_coef_names

stopifnot(identical(re_names2, expected_re))

stopifnot(identical(nrow(fit2$fixef[[re_names2[1L]]]), N_ING))



.validate_manual_block2_fixef(fit2, label = "ING estimated vcov", ing = TRUE)



.validate_lmerb_re(fit2, label = "ING estimated vcov")



cat("\ntest_lmerb_mer_re_validation: OK\n")



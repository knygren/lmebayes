test_that("lmerb reuses Prior_Setup_lmebayes()'s glmmTMB fit instead of re-fitting it", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(
    form,
    data            = dat,
    pwt             = 0.01,
    pwt_measurement = 0.1,
    dispformula     = ~Subject
  )
  disp_pf <- dGamma_list(ps, warn_asymmetric = FALSE)
  expect_identical(attr(disp_pf, "dispersion_fit"), ps$dispersion_fit)

  set.seed(1L)
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = disp_pf,
    dispformula      = ~Subject,
    n                = 50L
  )

  expect_identical(fit$prior$dispersion_mode, "gamma_list")
  expect_s4_class(fit$lmer, "merMod")
  ## Reused ps$dispersion_fit, not a second glmmTMB fit of the same model.
  expect_identical(fit$dispersion_fit, ps$dispersion_fit)
})

test_that("dispformula validation rejects mismatched dispersion_ranef shapes", {
  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)

  expect_error(
    lmerb(
      form,
      data             = dat,
      pfamily_list     = pfamily_list(ps),
      dispersion_ranef = ps$dispersion_ranef,
      dispformula      = ~Subject,
      n                = 2L
    ),
    "dGamma_list"
  )

  skip_if_not_installed("glmmTMB")
  ps_grp <- Prior_Setup_lmebayes(
    form,
    data            = dat,
    pwt             = 0.01,
    pwt_measurement = 0.1,
    dispformula     = ~Subject
  )
  disp_pf <- dGamma_list(ps_grp, warn_asymmetric = FALSE)
  expect_error(
    lmerb(
      form,
      data             = dat,
      pfamily_list     = pfamily_list(ps_grp),
      dispersion_ranef = disp_pf,
      dispformula      = ~1,
      n                = 2L
    ),
    "dispformula"
  )
})

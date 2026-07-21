## sim_method wiring for lmerb()/glmerb(): the fixed-dispersion /
## known-variance-components route (all-dNormal pfamily_list + fixed scalar
## or fixed per-group dispersion_ranef) has both an exact-iid engine
## ("DEFAULT") and the two-block Gibbs engine ("TWO_BLOCK_GIBBS"); every
## other route only has Gibbs, so sim_method is a no-op there. Exactness of
## the iid engine itself is covered in lmebayesCore's
## test-sim-method-iid.R; this file only checks that lmerb()/glmerb() thread
## the argument through correctly and that print()/summary() report the
## engine that actually ran.

test_that("lmerb(): sim_method dispatches DEFAULT to iid and TWO_BLOCK_GIBBS to Gibbs", {
  skip_on_cran()

  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)

  set.seed(1L)
  fit_default <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 50L,
    progbar          = FALSE
  )
  expect_identical(fit_default$sim_method_used, "DEFAULT")
  expect_identical(fit_default$draw_engine, "rLMMNormal_joint_iid")
  expect_output(print(fit_default), "exact iid")
  expect_output(print(summary(fit_default)), "exact iid")

  set.seed(1L)
  fit_gibbs <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 50L,
    progbar          = FALSE,
    sim_method       = "TWO_BLOCK_GIBBS"
  )
  expect_identical(fit_gibbs$sim_method_used, "TWO_BLOCK_GIBBS")
  expect_false(identical(fit_gibbs$draw_engine, "rLMMNormal_joint_iid"))
  expect_output(print(fit_gibbs), "two-block Gibbs")
  expect_output(print(summary(fit_gibbs)), "two-block Gibbs")

  ## Both engines target the same exact posterior; fixef.mode (the exact
  ## ICM mean) must agree regardless of sim_method.
  expect_equal(fit_default$fixef.mode, fit_gibbs$fixef.mode)
})

test_that("lmerb(): sim_method is a no-op off the fixed+known-vcov route", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(big_word_club, package = "bayesrules", envir = environment())
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(dat, !is.na(score_ppvt))
  form <- score_ppvt ~ private_school + (1 | school_id)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
  pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

  set.seed(1L)
  fit_default <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pf,
    dispersion_ranef = ps$dispersion_ranef,
    n                = 20L,
    gap_tol          = NULL,
    progbar          = FALSE
  )
  expect_identical(fit_default$sim_method_used, "TWO_BLOCK_GIBBS")

  set.seed(1L)
  fit_explicit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pf,
    dispersion_ranef = ps$dispersion_ranef,
    n                = 20L,
    gap_tol          = NULL,
    progbar          = FALSE,
    sim_method       = "TWO_BLOCK_GIBBS"
  )
  expect_identical(fit_explicit$sim_method_used, "TWO_BLOCK_GIBBS")
  expect_identical(fit_default$draw_engine, fit_explicit$draw_engine)
  expect_output(print(fit_default), "two-block Gibbs")
})

test_that("lmerb()/glmerb() reject unknown sim_method values", {
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
      n                = 2L,
      sim_method       = "bogus"
    ),
    "sim_method"
  )

  expect_error(
    glmerb(
      form,
      data             = dat,
      family           = gaussian(),
      pfamily_list     = pfamily_list(ps),
      dispersion_ranef = ps$dispersion_ranef,
      n                = 2L,
      sim_method       = "bogus"
    ),
    "sim_method"
  )
})

test_that("glmerb(family = gaussian()): sim_method dispatches like lmerb()", {
  skip_on_cran()

  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)

  set.seed(1L)
  fit_default <- glmerb(
    form,
    data             = dat,
    family           = gaussian(),
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 50L,
    progbar          = FALSE
  )
  expect_identical(fit_default$sim_method_used, "DEFAULT")
  expect_identical(fit_default$draw_engine, "rLMMNormal_joint_iid")
  expect_output(print(fit_default), "exact iid")

  set.seed(1L)
  fit_gibbs <- glmerb(
    form,
    data             = dat,
    family           = gaussian(),
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 50L,
    progbar          = FALSE,
    sim_method       = "TWO_BLOCK_GIBBS"
  )
  expect_identical(fit_gibbs$sim_method_used, "TWO_BLOCK_GIBBS")
  expect_output(print(fit_gibbs), "two-block Gibbs")
})

test_that("lmerb(): simulate = FALSE leaves sim_method_used NULL", {
  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)

  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    simulate         = FALSE
  )
  expect_null(fit$sim_method_used)
})

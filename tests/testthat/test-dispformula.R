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

test_that("lmerb accepts a fixed named numeric vector for dispersion_ranef (fixed_vector mode)", {
  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
  grp_levels <- levels(dat$Subject)
  set.seed(2L)
  disp_vec <- stats::setNames(
    as.numeric(ps$dispersion_ranef) * stats::runif(length(grp_levels), 0.5, 1.5),
    sample(grp_levels)
  )

  set.seed(1L)
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = disp_vec,
    dispformula      = ~Subject,
    n                = 20L,
    progbar          = FALSE
  )

  expect_identical(fit$prior$dispersion_mode, "fixed_vector")
  expect_null(fit$dispersion_fit)
  expect_named(fit$sigma2, grp_levels)
  expect_equal(as.numeric(fit$sigma2), as.numeric(disp_vec[grp_levels]))
  expect_identical(fit$sigma2, fit$sigma2.mean)

  summ <- summary_sigma2(fit)
  expect_s3_class(summ, "summary.sigma2.lmerb")
  expect_identical(summ$mode, "fixed_vector")
  expect_null(summ$overview)
  expect_equal(
    summ$prior[grp_levels, "Fixed sigma2"],
    as.numeric(disp_vec[grp_levels])
  )
})

test_that("dispersion_ranef fixed_vector validation rejects bad shapes and ~1", {
  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
  grp_levels <- levels(dat$Subject)

  ## Wrong length.
  expect_error(
    lmerb(
      form,
      data             = dat,
      pfamily_list     = pfamily_list(ps),
      dispersion_ranef = stats::setNames(rep(1, 3), grp_levels[1:3]),
      dispformula      = ~Subject,
      n                = 2L
    ),
    "length"
  )

  ## Mismatched names.
  bad_vec <- stats::setNames(rep(1, length(grp_levels)), grp_levels)
  names(bad_vec)[1] <- "not_a_subject"
  expect_error(
    lmerb(
      form,
      data             = dat,
      pfamily_list     = pfamily_list(ps),
      dispersion_ranef = bad_vec,
      dispformula      = ~Subject,
      n                = 2L
    ),
    "group levels"
  )

  ## dispformula = ~1 with a per-group vector is rejected.
  good_vec <- stats::setNames(rep(1, length(grp_levels)), grp_levels)
  expect_error(
    lmerb(
      form,
      data             = dat,
      pfamily_list     = pfamily_list(ps),
      dispersion_ranef = good_vec,
      dispformula      = ~1,
      n                = 2L
    ),
    "dispformula"
  )
})

test_that("glmerb(family = gaussian()) matches lmerb() for fixed_vector dispersion_ranef", {
  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
  grp_levels <- levels(dat$Subject)
  set.seed(3L)
  disp_vec <- stats::setNames(
    as.numeric(ps$dispersion_ranef) * stats::runif(length(grp_levels), 0.5, 1.5),
    grp_levels
  )

  set.seed(1L)
  fit_glm <- glmerb(
    form,
    data             = dat,
    family           = gaussian(),
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = disp_vec,
    dispformula      = ~Subject,
    n                = 20L,
    progbar          = FALSE
  )

  expect_identical(fit_glm$prior$dispersion_mode, "fixed_vector")
  expect_null(fit_glm$dispersion_fit)
  expect_named(fit_glm$sigma2, grp_levels)
  expect_equal(as.numeric(fit_glm$sigma2), as.numeric(disp_vec[grp_levels]))
})

test_that("summary.lmerb overview includes glmer reference and Pr(Prior_tail)", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)
  set.seed(1L)
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 50L
  )

  sm <- summary(fit)
  ov <- sm$fixef_overview
  pov <- sm$fixef_prior_overview
  expect_equal(rownames(pov), rownames(ov))
  expect_true(all(c("Prior Mean", "Prior.sd", "lmer", "lmer.se") %in% colnames(pov)))
  expect_false(any(c("lmer", "lmer.se") %in% colnames(ov)))
  expect_true("Pr(Prior_tail)" %in% colnames(ov))
  expect_true(all(is.finite(ov[["Pr(Prior_tail)"]])))

  ref <- sm$fixef[["(Intercept)"]]$coefficients1
  expect_true(all(c("lmer", "lmer.se") %in% colnames(ref)))
  expect_equal(
    pov["(Intercept)::(Intercept)", "lmer"],
    ref["(Intercept)", "lmer"]
  )

  pct <- sm$fixef_percentiles_overview
  expect_equal(rownames(pct), rownames(ov))
  expect_true(all(c("1.0%", "Median", "99.0%") %in% colnames(pct)))
  expect_equal(
    pct["(Intercept)::(Intercept)", "Median"],
    sm$fixef[["(Intercept)"]]$Percentiles["(Intercept)", "Median"]
  )

  ro <- sm$ranef_overview
  cand_col <- grep("^Cand/draw$", colnames(ro), value = TRUE)
  expect_length(cand_col, 0L)

  expect_true(is.null(sm$ranef.iters.mean))
})

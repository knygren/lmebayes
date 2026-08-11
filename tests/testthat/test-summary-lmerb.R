test_that("summary.lmerb overview includes glmer reference and Pr(Prior_tail)", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
  set.seed(1L)
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$group.dispersion,
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

test_that("ranef/coef/fixef mirror lme4 layout on lmerb fit", {
  skip_on_cran()

  data(sleepstudy, package = "lme4", envir = environment())
  dat <- sleepstudy
  dat$Subject <- factor(dat$Subject)
  form <- Reaction ~ Days + (Days || Subject)

  ps <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)
  set.seed(2L)
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$group.dispersion,
    n                = 40L
  )

  re_mode <- lme4::ranef(fit, type = "mode")
  expect_s3_class(re_mode, "ranef.lmerb")
  expect_equal(
    as.matrix(re_mode$Subject),
    fit$groupef.mode,
    tolerance = 1e-10
  )

  re_mean <- lme4::ranef(fit, type = "mean", condVar = TRUE)
  expect_true(!is.null(attr(re_mean, "postVar")))
  manual_mean <- tapply(
    seq_len(nrow(fit$groupef)),
    fit$groupef$Subject,
    function(idx) colMeans(fit$groupef[idx, fit$model_setup$groupef.names, drop = FALSE]),
    simplify = FALSE
  )
  manual_mat <- do.call(rbind, manual_mean[rownames(fit$groupef.mode)])
  expect_equal(as.matrix(re_mean$Subject), manual_mat, tolerance = 1e-10)

  long <- as.data.frame(re_mean)
  expect_true(all(c("grpvar", "term", "grp", "condval", "condsd") %in% names(long)))
  expect_true(all(is.finite(long$condsd)))

  expect_equal(unname(lme4::fixef(fit)), unname(lme4::fixef(fit$lmer)))

  cf <- coef(fit, type = "mode")
  expect_s3_class(cf, "coef.lmerb")
  mer_cf <- as.matrix(coef(fit$lmer)$Subject)
  expect_equal(as.matrix(cf$Subject), mer_cf, tolerance = 0.25)

  cf_mean <- coef(fit, type = "mean")
  expect_equal(as.matrix(cf_mean$Subject), manual_mat, tolerance = 1e-10)
})

# lmerb on bayesrules::big_word_club -- full cross-level moderation model
# (same as demo/Ex_12_lmerb_BigWordClub.R), plus a statistical check of the
# sampler against the exact Gaussian posterior.
#
# With dNormal pfamilies and fixed dispersions the joint posterior is
# Gaussian, so the ICM fixed point coef.mode is the *exact* posterior mean
# (= mode), and it is also the sampler's starting state.  The MCMC mean of
# the Block 2 draws must therefore agree with coef.mode up to Monte Carlo
# error: z = (draws mean - coef.mode) / (SD / sqrt(n)) is approximately
# N(0, 1) per coefficient (stored draws are near-independent by the
# m_convergence TV calibration).  |z| < 4 keeps the false-failure
# probability negligible across all hyperparameters while still catching a
# drifting or mis-centered sampler.

test_that("lmerb: simulated Block 2 means match the exact posterior mean", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(big_word_club, package = "bayesrules", envir = environment())
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

  form <- score_ppvt ~
    private_school + title1 + free_reduced_lunch +
    distracted_ppvt + distracted_a1 +
    free_reduced_lunch:distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)

  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)

  set.seed(1L)
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 1000L
  )

  expect_s3_class(fit, "lmerb")
  re_names <- fit$model_setup$re_coef_names
  expect_identical(re_names, c("(Intercept)", "distracted_ppvt", "distracted_a1"))
  n_draws <- nrow(fit$fixef[[re_names[1L]]])
  expect_identical(n_draws, 1000L)

  for (k in re_names) {
    draws <- fit$fixef[[k]]
    dm    <- fit$fixef.means[[k]]
    icm   <- fit$fixef.mode[[k]]
    sd_k  <- apply(draws, 2L, sd)
    expect_true(all(is.finite(draws)))
    expect_true(all(sd_k > 0))

    z <- (dm - icm) / (sd_k / sqrt(n_draws))
    expect_true(
      all(abs(z) < 4),
      info = sprintf(
        "[%s] |z| too large: %s", k,
        paste(sprintf("%s = %.2f", names(z), z), collapse = ", ")
      )
    )
  }

  expect_no_error(capture.output(print(summary(fit))))
})

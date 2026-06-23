# glmerb Poisson smoke test on bayesrules::airbnb (full cross-level model;
# same as demo/Ex_13_glmerb_Airbnb.R and data-raw/test_glmerb_airbnb.R).
#
# Non-Gaussian glmerb runs a pilot stage from the ICM mode, then main-stage
# sampling from the pilot mean.  For Poisson GLMMs the posterior is skewed, so
# main-stage means need not match the ICM mode; we instead check that the
# two-stage workflow completes, pilot metadata is recorded, the hyperparameter
# design is nontrivial, and group-level random-effect summaries track
# coef(glmer) (ordering / scramble diagnostics).

test_that("glmerb: Poisson pilot workflow and group ordering on full airbnb", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(airbnb, package = "bayesrules", envir = environment())
  dat <- airbnb
  dat$rating_c    <- dat$rating - mean(dat$rating)
  dat$log_price_c <- scale(log(dat$price + 1))[, 1]
  dat$walk_c      <- dat$walk_score - mean(dat$walk_score)
  dat$transit_c   <- dat$transit_score - mean(dat$transit_score)
  dat <- dat[complete.cases(dat[, c(
    "reviews", "rating", "rating_c", "price", "log_price_c",
    "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
  )]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))

  form <- reviews ~
    walk_c + transit_c +
    rating_c + log_price_c +
    walk_c:rating_c + transit_c:log_price_c +
    (1 + rating_c + log_price_c || neighborhood)

  design <- model_setup(form, data = dat, family = poisson())
  expect_true(design$rank_ok)
  expect_identical(design$re_coef_names, c("(Intercept)", "rating_c", "log_price_c"))
  expect_gte(ncol(design$X_hyper[["(Intercept)"]]), 3L)
  expect_gte(ncol(design$X_hyper[["rating_c"]]), 2L)
  expect_gte(ncol(design$X_hyper[["log_price_c"]]), 2L)

  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)

  set.seed(42L)
  fit <- glmerb(
    form,
    data         = dat,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    n            = 1000L,
    mode_gap_max = 1.0
  )

  expect_s3_class(fit, "glmerb")
  re_names <- fit$model_setup$re_coef_names
  expect_identical(re_names, c("(Intercept)", "rating_c", "log_price_c"))
  n_draws <- nrow(fit$fixef[[re_names[1L]]])
  expect_identical(n_draws, 1000L)
  expect_false(is.null(fit$fixef.init))
  expect_type(fit$convergence, "list")
  expect_true(is.finite(fit$convergence$m_convergence))
  expect_gt(fit$pilot_chisq$n_pilot, 0L)
  expect_identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot)
  expect_identical(fit$convergence$n_pilot_source, "cost")
  expect_true(is.finite(fit$convergence$m_convergence_pilot))
  expect_identical(fit$convergence$mode_gap_max, 1.0)
  expect_true(is.finite(fit$pilot_chisq$p_value))
  expect_true(is.finite(fit$pilot_chisq$Q))

  grp_col  <- fit$model_setup$group_name
  grp_levs <- rownames(coef(fit$glmer)[[grp_col]])
  J        <- length(grp_levs)
  icm_b    <- fit$ranef.mode
  expect_true(setequal(rownames(icm_b), grp_levs))
  expect_true(setequal(colnames(fit$fixef.mu), grp_levs))

  re_draws_mean <- tapply(
    seq_len(nrow(fit$coefficients)),
    fit$coefficients[[grp_col]],
    function(idx) colMeans(fit$coefficients[idx, re_names, drop = FALSE]),
    simplify = FALSE
  )

  mcmc_mat <- do.call(rbind, lapply(grp_levs, function(l) {
    re_draws_mean[[as.character(l)]]
  }))
  rownames(mcmc_mat) <- grp_levs
  icm_mat <- icm_b[grp_levs, re_names, drop = FALSE]
  glmer_mat <- as.matrix(coef(fit$glmer)[[grp_col]][grp_levs, re_names, drop = FALSE])

  scl <- apply(icm_mat, 2L, sd)
  scl[scl < 1e-8] <- 1
  A <- sweep(mcmc_mat, 2L, scl, "/")
  B <- sweep(icm_mat, 2L, scl, "/")
  D <- outer(rowSums(A^2), rep(1, J)) +
    outer(rep(1, J), rowSums(B^2)) - 2 * A %*% t(B)
  nearest <- apply(D, 1L, which.min)
  n_match <- sum(nearest == seq_len(J))
  expect_gte(n_match, ceiling(0.9 * J))

  for (k in re_names) {
    expect_gt(cor(mcmc_mat[, k], glmer_mat[, k]), 0.8)
    expect_gt(cor(mcmc_mat[, k], icm_mat[, k]), 0.9)
  }

  expect_no_error(capture.output(print(summary(fit))))
})

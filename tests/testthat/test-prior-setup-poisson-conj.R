test_that("Prior_Setup conj_poisson matches pwt / n_prior conjugate identities", {
  n <- 9L
  y <- c(rep(1, 3), rep(0, 6))
  df <- data.frame(y = y)

  pwt <- 0.4
  ps <- Prior_Setup(
    y ~ 1,
    family = poisson(link = "identity"),
    data = df,
    weights = rep(1, n),
    pwt = pwt
  )

  expect_false(is.null(ps$conj_poisson))

  ybar <- mean(y)
  n_eff <- n
  np <- (pwt / (1 - pwt)) * n_eff

  expect_equal(ps$PriorSettings$n_prior, np)
  expect_equal(as.numeric(ps$conj_poisson$shape), np * ybar)
  expect_equal(as.numeric(ps$conj_poisson$rate), np)
  expect_equal(as.numeric(ps$conj_poisson$weighted_mean_rate), ybar)
  expect_equal(as.numeric(ps$conj_poisson$n_prior_eff), np)

  S <- sum(y)
  gam_shape <- ps$conj_poisson$shape + S
  gam_rate <- ps$conj_poisson$rate + sum(rep(1, n))
  expect_equal(gam_shape / gam_rate, ybar)
})

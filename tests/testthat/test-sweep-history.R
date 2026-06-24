# Sweep history is stored on glmerb fits and printed on demand (not during sampling).

test_that("glmerb stores sweep_history and print(sweep_history=TRUE) works", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(airbnb, package = "bayesrules", envir = environment())
  dat <- airbnb
  dat$rating_c <- dat$rating - mean(dat$rating)
  dat <- dat[complete.cases(dat[, c("reviews", "rating_c", "neighborhood")]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))

  form <- reviews ~ rating_c + (1 + rating_c || neighborhood)
  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)

  set.seed(1L)
  fit <- glmerb(
    form,
    data = dat,
    family = poisson(),
    pfamily_list = pfamily_list(ps),
    n = 50L,
    m_convergence = 5L,
    m_convergence_pilot = 3L,
    progbar = FALSE
  )

  expect_s3_class(fit$sweep_history$main, "two_block_sweep_history")
  expect_s3_class(fit$sweep_history$pilot, "two_block_sweep_history")
  expect_equal(fit$sweep_history$main$stage, "main")
  expect_equal(fit$sweep_history$pilot$stage, "pilot")

  out <- capture.output(
    print(fit, sweep_history = TRUE, sweep_history_stage = "main", max_sweeps = 2L)
  )
  expect_true(any(grepl("fixef by sweep", out, fixed = TRUE)))
})

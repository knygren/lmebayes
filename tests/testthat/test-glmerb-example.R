# glmerb (Poisson) on bayesrules::airbnb_small -- same model as the ?glmerb
# example, with an overall (multivariate) centering test.
#
# The sampler uses a pilot stage of independent chains from the ICM mode to
# estimate the posterior mean (coef.pilot.mean), then runs the main n draws
# from that center.
#
# We test centering overall (not coefficient-by-coefficient):
#   H0_pilot: E[beta | y] = coef.pilot.mean
#   H0_mode : E[beta | y] = coef.mode
# using Wald statistics with covariance(cbar) = cov(draws)/n.
#
# Target behaviour (p-value only):
#   - ideally, not significantly different from pilot mean (p > 0.05),
#   - if not, at least pilot null should be less rejected than mode null.

.print_glmerb_test_table <- function(x, digits = 6L) {
  op <- options(width = max(300L, getOption("width")))
  on.exit(options(op), add = TRUE)
  paste(
    capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE)),
    collapse = "\n"
  )
}

.print_glmerb_test_summary <- function(fit) {
  op <- options(width = max(300L, getOption("width")))
  on.exit(options(op), add = TRUE)
  paste(capture.output(print(summary(fit))), collapse = "\n")
}

test_that("glmerb: overall posterior mean is closer to pilot mean than mode", {
  skip_on_cran()
  skip_if_not_installed("bayesrules")

  data(airbnb_small, package = "bayesrules", envir = environment())
  dat <- airbnb_small
  dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
  dat$walk_c   <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
  dat <- dat[complete.cases(dat[, c("reviews", "rating_c", "walk_c",
                                    "neighborhood")]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))

  form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(),
                             pwt = 0.01)

  fit <- glmerb(
    form,
    data         = dat,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    n            = 2000L,
    n_pilot      = 2000L,
    seed         = 42L
  )

  expect_s3_class(fit, "glmerb")
  re_names <- fit$model_setup$re_coef_names
  expect_identical(re_names, c("(Intercept)", "rating_c"))

  n_draws <- nrow(fit$fixef_draws[[re_names[1L]]])
  expect_identical(n_draws, 2000L)

  # pilot mean and convergence metadata must have been computed
  expect_false(is.null(fit$coef.pilot.mean))
  expect_true(is.list(fit$convergence))
  expect_true(is.finite(fit$convergence$m_convergence))
  expect_identical(
    fit$convergence$m_convergence_pilot,
    fit$convergence$m_convergence
  )
  expect_true(is.list(fit$pilot_mode_test))
  expect_true(is.finite(fit$pilot_mode_test$p_value))
  expect_identical(fit$pilot_mode_test$n_pilot, 2000L)

  X <- do.call(cbind, lapply(re_names, function(k) fit$fixef_draws[[k]]))
  cn <- unlist(lapply(re_names, function(k) {
    paste0(k, "::", colnames(fit$fixef_draws[[k]]))
  }))
  colnames(X) <- cn

  expect_true(all(is.finite(X)))

  beta_bar <- colMeans(X)
  theta_pilot <- unlist(lapply(re_names, function(k) fit$coef.pilot.mean[[k]]))
  theta_mode  <- unlist(lapply(re_names, function(k) fit$coef.mode[[k]]))
  names(theta_pilot) <- cn
  names(theta_mode)  <- cn

  # Test-only diagnostic table (not user-facing):
  # side-by-side mode, pilot mean, and posterior mean from main draws.
  center_tab <- data.frame(
    parameter = cn,
    mode = unname(theta_mode),
    pilot_mean = unname(theta_pilot),
    main_mean = unname(beta_bar),
    stringsAsFactors = FALSE
  )
  rownames(center_tab) <- NULL
  diff_tab <- data.frame(
    parameter = cn,
    pilot_minus_mode = unname(theta_pilot - theta_mode),
    main_minus_pilot = unname(beta_bar - theta_pilot),
    main_minus_mode = unname(beta_bar - theta_mode),
    stringsAsFactors = FALSE
  )
  rownames(diff_tab) <- NULL

  sd_main <- apply(X, 2L, stats::sd)
  se_main <- sd_main / sqrt(n_draws)
  z_main_vs_pilot <- unname((beta_bar - theta_pilot) / se_main)
  z_main_vs_mode  <- unname((beta_bar - theta_mode) / se_main)
  uni_tab <- data.frame(
    parameter = cn,
    z_vs_pilot = z_main_vs_pilot,
    p_vs_pilot = 2 * stats::pnorm(abs(z_main_vs_pilot), lower.tail = FALSE),
    z_vs_mode = z_main_vs_mode,
    p_vs_mode = 2 * stats::pnorm(abs(z_main_vs_mode), lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
  rownames(uni_tab) <- NULL

  message("glmerb centers table (test-only):")
  message(.print_glmerb_test_table(center_tab))
  message("glmerb differences table (test-only):")
  message(.print_glmerb_test_table(diff_tab))
  message("glmerb univariate z/p table (test-only):")
  message(.print_glmerb_test_table(uni_tab))

  n_tot <- nrow(X)
  p_tot <- ncol(X)
  S <- stats::cov(X)
  V <- S / n_tot
  V_inv <- solve(V)

  d_pilot <- beta_bar - theta_pilot
  d_mode  <- beta_bar - theta_mode

  Q_pilot <- as.numeric(t(d_pilot) %*% V_inv %*% d_pilot)
  Q_mode  <- as.numeric(t(d_mode)  %*% V_inv %*% d_mode)
  p_pilot <- stats::pchisq(Q_pilot, df = p_tot, lower.tail = FALSE)
  p_mode  <- stats::pchisq(Q_mode,  df = p_tot, lower.tail = FALSE)

  message(sprintf(
    "glmerb centering test: m_convergence=%d, p(mean=pilot)=%.4g, p(mean=mode)=%.4g",
    fit$convergence$m_convergence, p_pilot, p_mode
  ))

  expect_true(
    is.finite(p_pilot) && is.finite(p_mode),
    info = "Overall p-values must be finite."
  )
  expect_true(
    (p_pilot > 0.05) || (p_pilot >= p_mode),
    info = sprintf(
      "Either p(mean=pilot) should exceed 0.05 or be less rejected than mode; got p_pilot=%.4g, p_mode=%.4g",
      p_pilot, p_mode
    )
  )

  out_print <- capture.output(print(fit))
  expect_false(any(grepl("m_convergence", out_print, fixed = TRUE)))
  expect_false(any(grepl("pilot_vs_mode", out_print, fixed = TRUE)))
  expect_false(any(grepl("pilot.mean", out_print, fixed = TRUE)))

  message("glmerb summary(fit) (test-only):")
  out_sum_text <- .print_glmerb_test_summary(fit)
  message(out_sum_text)
  out_sum <- strsplit(out_sum_text, "\n", fixed = TRUE)[[1L]]
  expect_s3_class(summary(fit), "summary.lmerb")
  expect_true(any(grepl("Block 2", out_sum, fixed = TRUE)))
  expect_true(any(grepl("coef.mode", out_sum, fixed = TRUE)))
  expect_false(any(grepl("m_convergence", out_sum, fixed = TRUE)))
  expect_false(any(grepl("Pilot vs mode", out_sum, fixed = TRUE)))
})

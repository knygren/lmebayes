# glmerb (Poisson) on a subset of bayesrules::airbnb_small, with an overall
# (multivariate) centering test.
#
# WHY A SUBSET:
#   The full airbnb_small has 17 neighborhood levels (J=17).  By a CLT-like
#   argument over J groups, the marginal posterior of the fixed-effect
#   hyperparameter gamma is approximately normal for large J, making the mode
#   nearly equal to the mean and rendering the pilot-vs-mode comparison
#   insensitive.  The commonly cited threshold is J ~ 30; at J=17 we are
#   already close to the normal regime.
#
#   To make the mode-vs-mean gap detectable we keep only the smallest
#   neighborhoods (5-20 observations, J=6, 72 rows total).
#   Two criteria amplify Poisson non-normality:
#     1. Few groups (J small) → posterior of gamma far from normal.
#     2. Small counts per observation → individual Poisson posteriors skewed;
#        ~43% of the review counts in these groups are single-digit.
#   Larger groups (>20 obs) and the very large ones (Logan Square n=330,
#   Rogers Park n=123) are excluded because for large n_j the Poisson
#   likelihood is well-approximated by a Gaussian, reducing skewness.
#   For Poisson (concave h): ICM mode < E[gamma|y] < gamma* (Banach fixed pt).
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

  # Subset to the smallest neighborhoods (<= 20 obs) to maximise Poisson
  # non-normality: few groups (small J) AND small per-observation counts.
  # Drop singletons and pairs (n=1,2) as they are too sparse to fit random
  # slopes reliably; keep groups with 5-20 observations.
  grp_counts  <- table(dat$neighborhood)
  keep_groups <- names(grp_counts[grp_counts >= 5L & grp_counts <= 20L])
  dat <- dat[dat$neighborhood %in% keep_groups, ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))
  # Sanity: expect exactly 6 groups (Montclare 5, O'Hare 5, North Park 7,
  # Lincoln Square 15, Jefferson Park 20, Portage Park 20) after subsetting.
  # 72 rows total, ~44% of review counts are single-digit.
  stopifnot(nlevels(dat$neighborhood) == 6L)

  form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(),
                             pwt = 0.01)

  fit <- glmerb(
    form,
    data         = dat,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    n            = 10000L,
    gap_tol      = 0.0196,  # default; gives n_pilot = 10000
    mode_gap_max = 1.0,     # default; gives m_convergence_pilot = 9 for p=3, lambda*~0.596
    seed         = 42L
  )

  expect_s3_class(fit, "glmerb")
  re_names <- fit$model_setup$re_coef_names
  expect_identical(re_names, c("(Intercept)", "rating_c"))

  n_draws <- nrow(fit$fixef_draws[[re_names[1L]]])
  expect_identical(n_draws, 10000L)

  # pilot mean and convergence metadata must have been computed
  expect_false(is.null(fit$coef.pilot.mean))
  expect_true(is.list(fit$convergence))
  expect_true(is.finite(fit$convergence$m_convergence))
  # gap_tol = 0.0196 => n_pilot = ceiling((qnorm(0.975)/0.0196)^2) = 10000
  expect_identical(fit$pilot_mode_test$n_pilot, 10000L)
  # mode_gap_max = 1, p = 3, lambda* ~ 0.596, tv_tol = 0.01:
  #   D_max = sqrt(3), c_tol = qnorm(0.505)/sqrt(2)*2*sqrt(2) ~ 0.0251
  #   l = ceil(log(1.732/0.0251)/log(1/0.596)) = ceil(8.19) = 9
  expect_identical(fit$convergence$m_convergence_pilot, 9L)
  expect_identical(fit$convergence$mode_gap_max, 1.0)
  expect_true(is.list(fit$pilot_mode_test))
  expect_true(is.finite(fit$pilot_mode_test$p_value))

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
  # "m_convergence" standalone (not part of "m_convergence_pilot") must not appear
  expect_false(any(grepl("m_convergence(?!_)", out_print, perl = TRUE)))
  expect_false(any(grepl("pilot_vs_mode", out_print, fixed = TRUE)))
  expect_false(any(grepl("pilot.mean", out_print, fixed = TRUE)))

  message("glmerb coef.means by component (test-only):")
  message(paste(capture.output(print_coef_means(fit)), collapse = "\n"))

  message("glmerb print(fit) (test-only):")
  message(paste(capture.output(print(fit)), collapse = "\n"))

  message("glmerb summary(fit) (test-only):")
  out_sum_text <- .print_glmerb_test_summary(fit)
  message(out_sum_text)
  out_sum <- strsplit(out_sum_text, "\n", fixed = TRUE)[[1L]]
  expect_s3_class(summary(fit), "summary.lmerb")
  expect_true(any(grepl("Block 2", out_sum, fixed = TRUE)))
  expect_true(any(grepl("coef.mode", out_sum, fixed = TRUE)))
  # "m_convergence" standalone (not part of "m_convergence_pilot") must not appear
  expect_false(any(grepl("m_convergence(?!_)", out_sum, perl = TRUE)))
  expect_false(any(grepl("Pilot vs mode", out_sum, fixed = TRUE)))
})

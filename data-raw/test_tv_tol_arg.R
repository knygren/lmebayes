# Regression test: tv_tol / m_convergence arguments on lmerb() and glmerb().
#
# lmerb / glmerb(gaussian): m_min is derived from the Theorem 3 TV bound
# (Nygren 2020) via two_block_rate() + two_block_l_for_tv(), + 1L for the
# half-step lag of the stored b draw.  glmerb(non-gaussian): same machinery
# on the local-Gaussian approximation at the ICM posterior mode
# (two_block_mode_weights), giving a LOWER BOUND m_min.  m_convergence
# (optional) overrides upward; values below m_min are raised with a warning.
#
#   Rscript data-raw/test_tv_tol_arg.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

data(big_word_club, package = "bayesrules")
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

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)

## --- expected m_convergence from the rate/bound machinery directly ---------
design <- model_setup(form_lmer, data = dat)
re_names <- design$re_coef_names
block1_prior <- lmebayes:::.lmebayes_block1_prior_list(ps)
block2_prior_list <- stats::setNames(
  lapply(re_names, function(k) {
    pl_k <- ps$prior_list[[k]]
    list(mu = pl_k$mu_fixef, Sigma = pl_k$Sigma_fixef,
         dispersion = pl_k$dispersion_fixef)
  }),
  re_names
)
rate <- glmbayesCore::two_block_rate(
  x = design$Z, block = design$groups, x_hyper = design$X_hyper,
  prior_list_block1 = block1_prior, prior_list_block2 = block2_prior_list,
  family = gaussian(), group_levels = levels(design$groups)
)
m_default <- glmbayesCore::two_block_l_for_tv(rate, 0.01) + 1L
m_strict  <- glmbayesCore::two_block_l_for_tv(rate, 1e-5) + 1L
stopifnot(m_strict > m_default)
cat(sprintf("expected m_convergence: tv_tol 0.01 -> %d, 1e-5 -> %d\n\n",
            m_default, m_strict))

## --- 1. lmerb default tv_tol = 0.01 ----------------------------------------
out1 <- capture.output(
  fit1 <- lmerb(form_lmer, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$dispersion_ranef,
                n = 25L, seed = 1L)
)
line1 <- grep("using m_convergence", out1, value = TRUE)
stopifnot(length(line1) == 1L)
stopifnot(grepl(sprintf("m_convergence = %d", m_default), line1))
stopifnot(inherits(fit1, "lmerb"))
stopifnot(identical(fit1$convergence$m_min, m_default))
stopifnot(identical(fit1$convergence$m_convergence, m_default))
stopifnot(identical(fit1$convergence$method, "exact"))
cat("1. lmerb default tv_tol:", line1, "\n")

## --- 2. lmerb stricter tv_tol ----------------------------------------------
out2 <- capture.output(
  fit2 <- lmerb(form_lmer, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$dispersion_ranef,
                n = 25L, tv_tol = 1e-5, seed = 1L)
)
line2 <- grep("using m_convergence", out2, value = TRUE)
stopifnot(grepl(sprintf("m_convergence = %d", m_strict), line2))
cat("2. lmerb tv_tol = 1e-5:", line2, "\n")

## --- 2b. m_convergence override: above m_min honored, below raised ----------
out2b <- capture.output(
  fit2b <- lmerb(form_lmer, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$dispersion_ranef,
                 n = 25L, m_convergence = 50L, seed = 1L)
)
stopifnot(identical(fit2b$convergence$m_convergence, 50L))
stopifnot(identical(fit2b$convergence$m_min, m_default))

warns <- character(0)
out2c <- withCallingHandlers(
  capture.output(
    fit2c <- lmerb(form_lmer, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$dispersion_ranef,
                   n = 25L, m_convergence = 2L, seed = 1L)
  ),
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(any(grepl("below the derived minimum", warns)))
stopifnot(identical(fit2c$convergence$m_convergence, m_default))
cat("2b. m_convergence override: 50 honored; 2 raised to m_min with warning\n")

## --- 3. invalid tv_tol errors ----------------------------------------------
for (bad in list(0, 1, -0.1, c(0.1, 0.2), "a", NA_real_)) {
  res <- tryCatch(
    lmerb(form_lmer, data = dat, pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$dispersion_ranef,
          n = 5L, tv_tol = bad),
    error = function(e) e
  )
  stopifnot(inherits(res, "error"), grepl("tv_tol", conditionMessage(res)))
}
cat("3. invalid tv_tol values rejected\n")

## --- 4. glmerb gaussian: same calibration as lmerb --------------------------
out4 <- capture.output(
  fit4 <- glmerb(form_lmer, data = dat, family = gaussian(),
                 pfamily_list = pfamily_list(ps),
                dispersion_ranef = ps$dispersion_ranef, n = 25L, seed = 1L)
)
line4 <- grep("using m_convergence", out4, value = TRUE)
stopifnot(length(line4) == 1L)
stopifnot(grepl(sprintf("m_convergence = %d", m_default), line4))
stopifnot(any(grepl("exact \\(Gaussian posterior\\)", out4)))
stopifnot(inherits(fit4, "glmerb"))
stopifnot(identical(fit4$convergence$method, "exact"))
cat("4. glmerb gaussian:", line4, "\n")

## --- 5. glmerb poisson: approximate local-Gaussian calibration --------------
data("airbnb_small", package = "bayesrules")
ab <- airbnb_small
ab$rating_c <- ab$rating - mean(ab$rating)
ab <- ab[complete.cases(ab[, c("reviews", "rating_c", "neighborhood")]), ]
form_pois <- reviews ~ rating_c + (1 + rating_c || neighborhood)
ps_pois <- Prior_Setup_lmebayes(form_pois, data = ab, family = poisson(),
                                pwt = 0.01)
out5 <- capture.output(
  fit5 <- glmerb(form_pois, data = ab, family = poisson(),
                 pfamily_list = pfamily_list(ps_pois), n = 10L,
                 tv_tol = 0.001, seed = 1L)
)
stopifnot(any(grepl("approximate \\(local-Gaussian at mode, poisson\\)", out5)))
stopifnot(inherits(fit5, "glmerb"))
stopifnot(identical(fit5$convergence$method, "local_gaussian_mode"))
stopifnot(fit5$convergence$m_min >= 1L)
stopifnot(identical(fit5$convergence$m_convergence, fit5$convergence$m_min))
stopifnot(fit5$convergence$lambda_star >= 0, fit5$convergence$lambda_star < 1)
cat(sprintf(
  "5. glmerb poisson: approximate calibration, lambda* = %.4f, m_min = %d\n",
  fit5$convergence$lambda_star, fit5$convergence$m_min
))

## --- 6. glmerb poisson: override doubles the lower bound --------------------
m_min5 <- fit5$convergence$m_min
out6 <- capture.output(
  fit6 <- glmerb(form_pois, data = ab, family = poisson(),
                 pfamily_list = pfamily_list(ps_pois), n = 10L,
                 tv_tol = 0.001, m_convergence = 2L * m_min5, seed = 1L)
)
stopifnot(identical(fit6$convergence$m_convergence, 2L * m_min5))
stopifnot(identical(fit6$convergence$m_min, m_min5))
cat(sprintf("6. glmerb poisson: override 2 * m_min = %d honored\n",
            2L * m_min5))

cat("\ntest_tv_tol_arg.R: all checks passed\n")

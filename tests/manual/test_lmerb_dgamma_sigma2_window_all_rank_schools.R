# Manual check: dGamma sigma^2 truncation window vs number of full-rank schools.
#
# Same RE structure as test_lmerb_dgamma_small5_validation.R (p_re = 2:
# Intercept + distracted_ppvt). Compares the 5-school subset to all algebraically
# full-rank schools on big_word_club for that model.
#
# With fixed pwt_measurement, n_prior grows with total n, so the 98% prior
# window on sigma^2 should tighten as k (and n) increase.
#
# Requires dev glmbayesCore for pwt_measurement.
#   cd .../lmebayes
#   Rscript tests/manual/test_lmerb_dgamma_sigma2_window_all_rank_schools.R

source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

PWT_MEASUREMENT <- 0.49

.print_sigma2_window <- function(ps, label) {
  m <- ps$ing_prior_measurement
  ratio <- m$disp_upper / m$disp_lower
  n_g <- length(ps$design$y)
  k   <- length(ps$design$re_rank)
  stopifnot(identical(k, sum(ps$design$re_rank)))
  cat(sprintf(
    paste0(
      "\n=== %s ===\n",
      "  schools (k)           : %d\n",
      "  observations (n)      : %d\n",
      "  mean obs/school       : %.2f\n",
      "  pwt_measurement       : %.2g\n",
      "  n_prior_measurement   : %.4g\n",
      "  shape / rate          : %.4g / %.4g\n",
      "  sigma2_hat (REML)     : %.4g\n",
      "  window [lower, upper] : [%.4g, %.4g]\n",
      "  ratio upper/lower     : %.4f\n",
      "  sigma2_hat / lower    : %.4f\n",
      "  upper / sigma2_hat    : %.4f\n\n"
    ),
    label,
    k,
    n_g,
    n_g / k,
    ps$pwt_measurement,
    ps$n_prior_measurement,
    m$shape,
    m$rate,
    m$sigma2_hat,
    m$disp_lower,
    m$disp_upper,
    ratio,
    m$sigma2_hat / m$disp_lower,
    m$disp_upper / m$sigma2_hat
  ))
  invisible(list(ratio = ratio, m = m, k = k, n = n_g))
}

cat("\n--- Reference: 5-school subset (same as small5 smoke test) ---\n")
fx5 <- .prepare_small5_lmerb_manual(n_schools = 5L)
ps5 <- Prior_Setup_lmebayes(
  fx5$form,
  data            = fx5$dat,
  pwt             = 0.01,
  pwt_measurement = PWT_MEASUREMENT
)
w5 <- .print_sigma2_window(ps5, "5 full-rank schools (small5 subset)")

cat("\n--- All full-rank schools (p_re = 2 model) ---\n")
fx_all <- .prepare_small5_all_full_rank_manual()
ps_all <- Prior_Setup_lmebayes(
  fx_all$form,
  data            = fx_all$dat,
  pwt             = 0.01,
  pwt_measurement = PWT_MEASUREMENT
)
w_all <- .print_sigma2_window(ps_all, "All full-rank schools")

cat(sprintf(
  paste0(
    "Summary: window ratio all-rank / 5-school = %.4f / %.4f = %.4f\n",
    "(values < 1 mean the all-rank window is tighter)\n"
  ),
  w_all$ratio,
  w5$ratio,
  w_all$ratio / w5$ratio
))

stopifnot(w_all$ratio < w5$ratio)
stopifnot(w_all$n > w5$n)
stopifnot(w_all$k > w5$k)

cat("\ntest_lmerb_dgamma_sigma2_window_all_rank_schools: OK\n")

library(lmebayes)
library(bayesrules)

data(book_banning, package = "bayesrules", envir = environment())
dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)

form_glmerb <- removed_i ~ violent_i + (1 + violent_i || state)
ps <- Prior_Setup_lmebayes(form_glmerb, data = dat, family = binomial(), pwt = 0.01)

set.seed(42L)
fit <- glmerb(
  form_glmerb,
  data         = dat,
  family       = binomial(),
  pfamily_list = pfamily_list(ps),
  n            = 500L,
  mode_gap_max = 1.0,
  progbar      = FALSE,
  verbose      = TRUE
)

mu <- unlist(lapply(fit$fixef_draws, colMeans))
mode <- unlist(fit$coef.mode)

cat("\n=== C++ path checks ===\n")
cat("draw_engine:", fit$convergence$draw_engine, "\n")
cat("m_convergence_pilot:", fit$convergence$m_convergence_pilot, "\n")
cat("m_convergence:", fit$convergence$m_convergence, "\n")
cat("n_pilot:", fit$convergence$n_pilot, "\n")
cat("fixef means:", paste(round(mu, 4), collapse = ", "), "\n")
cat("ICM mode:    ", paste(round(mode, 4), collapse = ", "), "\n")

stopifnot(identical(fit$convergence$draw_engine, "two_block_rNormal_reg_v5"))
stopifnot(identical(fit$convergence$m_convergence_pilot, 19L))
stopifnot(identical(fit$convergence$m_convergence, 11L))
stopifnot(abs(mu[1] - (-1.116)) < 0.08)
stopifnot(abs(mu[2] - 0.396) < 0.08)

cat("validate_cpp_book_banning: OK\n")

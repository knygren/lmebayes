library(lmebayes)
library(bayesrules)

data(book_banning, package = "bayesrules", envir = environment())
dat <- book_banning[, c("state", "removed", "violent")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)

form <- removed_i ~ violent_i + (1 + violent_i || state)
ps <- Prior_Setup_lmebayes(form, data = dat, family = binomial(), pwt = 0.01)

set.seed(42L)
fit <- glmerb(
  form,
  data         = dat,
  family       = binomial(),
  pfamily_list = pfamily_list(ps),
  n            = 20L,
  gap_tol      = 0.2,
  progbar      = FALSE
)

cat("draw_engine:", fit$convergence$draw_engine, "\n")
cat("n_pilot:", fit$convergence$n_pilot, "\n")

if (!identical(fit$convergence$draw_engine, "two_block_rNormal_reg_v5")) {
  stop("Expected two_block_rNormal_reg_v5 draw engine")
}

cat("smoke_glmerb_v5: OK\n")

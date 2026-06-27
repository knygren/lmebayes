## Ex_16 book-banning run (C++ engine); progbar=FALSE for log capture
suppressPackageStartupMessages({
  library(lmebayes)
  library(bayesrules)
  library(lme4)
})

data(book_banning, package = "bayesrules", envir = environment())

dat <- book_banning[, c(
  "state", "removed", "violent", "antifamily", "language"
)]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)

form_book <- removed ~ violent + antifamily + language + (1 | state)
form_glmerb <- removed_i ~ violent_i + (1 + violent_i || state)

ps <- Prior_Setup_lmebayes(form_glmerb, data = dat, family = binomial(), pwt = 0.01)

cat("draw_engine constant:", lmebayes:::.rglmerb_engine, "\n\n")

set.seed(42L)
fit <- glmerb(
  form_glmerb,
  data         = dat,
  family       = binomial(),
  pfamily_list = pfamily_list(ps),
  n            = 3000L,
  mode_gap_max = 1.0,
  progbar      = FALSE,
  verbose      = TRUE
)

cat("\n=== Post-fit checks ===\n")
cat("draw_engine:", fit$convergence$draw_engine, "\n")
cat("m_convergence_pilot:", fit$convergence$m_convergence_pilot, "\n")
cat("m_convergence:", fit$convergence$m_convergence, "\n")
cat("n_pilot:", fit$convergence$n_pilot, "\n")

cat("\nm_convergence used:", fit$convergence$m_convergence, "\n")
cat(sprintf(
  "Pilot vs mode (chi-squared): p = %.4g (n_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot
))

cat("\n--- Ch. 18 reference glmer (random intercept; all three reasons) ---\n")
fit_book <- lme4::glmer(form_book, data = dat, family = binomial())
print(lme4::fixef(fit_book))

lmebayes:::print_coef_means(fit)
print(fit)
summary(fit)

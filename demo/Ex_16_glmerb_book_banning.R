## Demo: glmerb() binomial workflow on bayesrules::book_banning
##
## Bayes Rules! Ch. 18.2--18.4 hierarchical logistic regression of book-removal
## outcomes on challenge reasons with state random effects.  The textbook uses
##
##   removed ~ violent + antifamily + language + (1 | state)
##
## lmebayes requires climber-/challenge-level predictors to appear as population
## mean slopes with matching random effects (see demo/Ex_13_glmerb_Airbnb.R).
## This demo fits a minimal binomial glmerb model with violent coded as 0/1
## (violent_i) and a random slope on violent_i, and compares the classical
## glmer fit for the Ch. 18 random-intercept formula.
##
##   demo("Ex_16_glmerb_book_banning", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

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

print(ps)

set.seed(42L)
fit <- glmerb(
  form_glmerb,
  data         = dat,
  family       = binomial(),
  pfamily_list = pfamily_list(ps),
  n            = 1000L,
  mode_gap_max = 1.0,
  progbar = TRUE
)

cat("m_convergence used:", fit$convergence$m_convergence, "\n")
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

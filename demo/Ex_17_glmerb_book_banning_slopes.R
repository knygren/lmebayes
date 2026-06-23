## Demo: glmerb() binomial workflow on bayesrules::book_banning (cross-level slopes)
##
## Extends demo/Ex_16_glmerb_book_banning.R: state-level median income moderates
## the random slope for violent challenges (same cross-level pattern as
## demo/Ex_13_glmerb_Airbnb.R and demo/Ex_12_lmerb_BigWordClub.R).
##
##   removed_i ~ violent_i + income_c + income_c:violent_i +
##     (1 + violent_i || state)
##
##   demo("Ex_17_glmerb_book_banning_slopes", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

data(book_banning, package = "bayesrules", envir = environment())

dat <- book_banning[, c(
  "state", "removed", "violent", "median_income"
)]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
dat$income_c <- as.numeric(scale(dat$median_income))

## Reason-level removal summary: demo("Ex_18_glmerb_book_banning_all_reasons")

form_glmerb <- removed_i ~ violent_i + income_c + income_c:violent_i +
  (1 + violent_i || state)

design <- model_setup(form_glmerb, data = dat, family = binomial())
cat("\n=== model_setup (X_hyper per RE component) ===\n\n")
for (k in design$re_coef_names) {
  cat(k, ":\n")
  print(colnames(design$X_hyper[[k]]))
}

ps <- Prior_Setup_lmebayes(form_glmerb, data = dat, family = binomial(), pwt = 0.01)

print(ps)

set.seed(42L)
fit <- glmerb(
  form_glmerb,
  data         = dat,
  family       = binomial(),
  pfamily_list = pfamily_list(ps),
  n            = 1000L,
  mode_gap_max = 1.0
)

cat("\nm_convergence used:", fit$convergence$m_convergence, "\n")
cat(sprintf(
  "Pilot vs mode (chi-squared): p = %.4g (n_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot
))

lmebayes:::print_coef_means(fit)
print(fit)
summary(fit)

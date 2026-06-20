## lmerb: school random effects on bayesrules::big_word_club
##
## Random intercept and distracted_ppvt random slope by school; the school's
## share of free/reduced-lunch students (free_reduced_lunch, constant within
## school) is a level-2 predictor of the school intercepts.
##
## The full workflow with factor-level diagnostics is preserved as a demo:
##   demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")

if (requireNamespace("bayesrules", quietly = TRUE)) {

  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(dat, !is.na(score_ppvt) & !is.na(invalid_ppvt) &
                  invalid_ppvt == 0L)
  dat <- dat[complete.cases(dat[, c("score_ppvt", "distracted_ppvt",
                                    "free_reduced_lunch", "school_id")]), ]

  form <- score_ppvt ~ free_reduced_lunch + distracted_ppvt +
    (1 + distracted_ppvt || school_id)

  ## Default hyperpriors calibrated from a reference lmer fit (weak prior).
  ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01)

  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pfamily_list(ps),
    dispersion_ranef = ps$dispersion_ranef,
    n                = 1000L,
    seed             = 1L
  )
  ## Block 2 posterior means alongside lmer MLE and ICM posterior mean.
  lmebayes:::print_coef_means(fit)
  print(fit)
  summary(fit)
}

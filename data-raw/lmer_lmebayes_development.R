## lmerb -- development script
##
## Draft entry point for a naive / empirical-Bayes lmerb() fit on the
## bayesrules::big_word_club school-based example (same model as
## Prior_Setup_lmebayes_development.R).
##
## Current lmerb() implementation: call model_setup() and return its
## output.  Later steps will wire in Prior_Setup_lmebayes(), Block 1
## (block_rNormalReg), and Block 2 (rNormal_reg per RE coefficient).

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules: install.packages('bayesrules')")
}
pkgload::load_all(export_all = FALSE)

## ===========================================================================
## Development run: big_word_club with free_reduced_lunch x distracted_a1
## ===========================================================================

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat$age_c     <- dat$age_months - mean(dat$age_months, na.rm = TRUE)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c("score_ppvt", "age_c", "distracted_a1",
                           "private_school", "title1", "free_reduced_lunch",
                           "school_id")])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_a1 + free_reduced_lunch:distracted_a1 +
  (1 + age_c + distracted_a1 || school_id)

ctrl_bobyqa <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

cat("=== lmerb (draft): model_setup only ===\n\n")
fit <- lmerb(
  form_lmer,
  data = dat,
  measurement_prior_list = NULL,
  control = ctrl_bobyqa
)
print(fit)

cat("\n=== Reference: lmer on same formula ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, control = ctrl_bobyqa)
print(summary(fit_lmer))

cat("\n=== Quick check: lmerb lmer fit matches standalone lmer ===\n")
stopifnot(all.equal(lme4::fixef(fit$lmer_fit), lme4::fixef(fit_lmer), check.attributes = FALSE))
stopifnot(all.equal(as.matrix(lme4::VarCorr(fit$lmer_fit)), as.matrix(lme4::VarCorr(fit_lmer)), tolerance = 1e-6))

cat("\nlmerb development script: OK\n")

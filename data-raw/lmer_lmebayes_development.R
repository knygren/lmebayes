## lmerb -- development script
##
## Block 1 smoke test on bayesrules::big_word_club (school-based example).
## Prior_Setup_lmebayes() must be called explicitly; lmerb() does not
## construct priors internally.

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules: install.packages('bayesrules')")
}
pkgload::load_all(export_all = FALSE)

## ===========================================================================
## Data and formula (same as Prior_Setup_lmebayes_development.R)
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

## ===========================================================================
## Step 1: explicit prior setup (required before lmerb)
## ===========================================================================

cat("=== Prior_Setup_lmebayes ===\n\n")
ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)

## ===========================================================================
## Step 2: lmerb Block 1 draw
## ===========================================================================

cat("\n=== lmerb: Block 1 random-effects draw ===\n\n")
fit <- lmerb(
  form_lmer,
  ps,
  n = 1L,
  data = dat,
  seed = 42L
)

cat(sprintf("  class: %s\n", paste(class(fit), collapse = ", ")))
cat(sprintf("  names: %s\n", paste(names(fit), collapse = ", ")))
cat(sprintf("  coefficients: %d rows x %d cols\n",
            nrow(fit$coefficients), ncol(fit$coefficients)))
cat("  columns:", paste(names(fit$coefficients), collapse = ", "), "\n")
cat("\n  first 5 rows:\n")
print(head(fit$coefficients, 5L))

## ===========================================================================
## Checks
## ===========================================================================

stopifnot(inherits(fit, "lmerb"))
stopifnot(identical(names(fit), c("model_setup", "lmer", "coefficients")))
stopifnot(inherits(fit$model_setup, "model_setup"))
stopifnot(inherits(fit$lmer, "lmerMod"))
stopifnot(nrow(fit$coefficients) == nlevels(ps$design$groups))
stopifnot(identical(
  names(fit$coefficients),
  c("draw", ps$design$group_name, ps$design$re_coef_names)
))

cat("\n=== Reference: lmer fixed effects (same formula) ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, control = ctrl_bobyqa)
print(lme4::fixef(fit_lmer))

cat("\nlmerb Block 1 development script: OK\n")

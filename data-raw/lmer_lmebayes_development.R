## lmerb -- development script
##
## Block 1 smoke test on bayesrules::big_word_club.
## Same formula and data as inst/examples/Ex_lmerb.R.

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
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c("score_ppvt", "distracted_a1", "distracted_ppvt",
                           "private_school", "title1", "free_reduced_lunch",
                           "school_id")])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_a1 + distracted_ppvt +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

ctrl_bobyqa <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## ===========================================================================
## Step 1: prior setup (required before lmerb)
## ===========================================================================

cat("=== Prior_Setup_lmebayes ===\n\n")
ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)

## ===========================================================================
## Step 2: lmerb — lmer on formula/data + Block 1 draws
## ===========================================================================

cat("\n=== lmerb: lmer fit + Block 1 random-effects draws ===\n\n")
fit <- lmerb(
  form_lmer,
  data = dat,
  measurement_prior_list = ps,
  n = 1L,
  control = ctrl_bobyqa,
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

design <- fit$model_setup
stopifnot(inherits(fit, "lmerb"))
stopifnot(identical(names(fit), c("model_setup", "lmer", "mu_all", "fixef",
                                   "fixef_mean", "fixef_draws", "fixef_draws_mean",
                                   "coefficients")))
stopifnot(all(dim(fit$mu_all) == c(length(design$re_coef_names), nlevels(design$groups))))
stopifnot(inherits(fit$model_setup, "model_setup"))
stopifnot(inherits(fit$lmer, "lmerMod"))
stopifnot(nrow(fit$coefficients) == nlevels(design$groups))
stopifnot(identical(
  names(fit$coefficients),
  c("draw", design$group_name, design$re_coef_names)
))

cat("\n=== Reference: standalone lmer fixed effects (same formula/data) ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, control = ctrl_bobyqa)
print(lme4::fixef(fit_lmer))
stopifnot(identical(lme4::fixef(fit$lmer), lme4::fixef(fit_lmer)))

cat("\nlmerb Block 1 development script: OK\n")

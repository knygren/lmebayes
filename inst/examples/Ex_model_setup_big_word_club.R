## model_setup() on bayesrules::big_word_club
##
## Same model as the full lmerb demo, demo/Ex_12_lmerb_BigWordClub.R (see
## also Prior_Setup_lmebayes / lmerb development scripts in data-raw/).
##
## Level 1 (students):
##   y ~ b0[j] + b_ppvt[j]*distracted_ppvt + b_a1[j]*distracted_a1
##
## Level 2 (schools):
##   b0[j]      ~ private_school + title1 + free_reduced_lunch + u0[j]
##   b_ppvt[j]  ~ 1 + u_ppvt[j]
##   b_a1[j]    ~ 1 + free_reduced_lunch + u_a1[j]
##                 (cross-level: free_reduced_lunch:distracted_a1 in formula)
##
## Each random slope has a matching fixed main effect (required for
## Prior_Setup_lmebayes() default calibration).

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_a1 + distracted_ppvt +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

ctrl_bobyqa <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## ---------------------------------------------------------------------------
## 1. lmer fit: raw output
## ---------------------------------------------------------------------------
cat("--- lmer fit ---\n")
fit <- lme4::lmer(form_lmer, data = dat, control = ctrl_bobyqa)
print(summary(fit))

cat("\n--- fixef(fit): population-level (gamma) estimates ---\n")
print(lme4::fixef(fit))

cat("\n--- coef(fit): per-group coefficients (fixef + ranef) ---\n")
print(coef(fit))

## ---------------------------------------------------------------------------
## 2. model_setup: structured view of the same model
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat, control = ctrl_bobyqa)
print(design)

## ---------------------------------------------------------------------------
## 3. Random effects b[j]: first 10 schools
## ---------------------------------------------------------------------------
cat("--- Random effects b[j]: first 10", design$group_name, "---\n")
re_df <- as.data.frame(lme4::ranef(design$lmer_fit)[[design$group_name]])
print(utils::head(re_df, 10))

## ---------------------------------------------------------------------------
## 4. Gamma estimates organised to match the Random Effects Model above
##
##    Mapping (intercept RE): X_hyper columns map directly to fixef() names.
##    Mapping (slope RE, hyper ~ 1): (Intercept) column -> fixef[slope_name].
## ---------------------------------------------------------------------------
cat("\n--- Random effects model (gamma estimates) ---\n")

fe         <- lme4::fixef(design$lmer_fit)
coef_df    <- coef(design$lmer_fit)[[design$group_name]]
coef_means <- colMeans(coef_df)
coef_vars  <- apply(coef_df, 2L, var)
coef_sds   <- sqrt(coef_vars)
w          <- max(nchar(design$re_coef_names))

for (nm in design$re_coef_names) {
  Xj    <- design$X_hyper[[nm]]
  other <- setdiff(colnames(Xj), "(Intercept)")
  hyper_rhs <- if (length(other) == 0L) "1" else paste(c("1", other), collapse = " + ")

  gamma <- setNames(
    vapply(colnames(Xj), function(col) {
      if (nm == "(Intercept)") {
        if (col %in% names(fe)) unname(fe[col]) else 0
      } else if (col == "(Intercept)") {
        if (nm %in% names(fe)) unname(fe[nm]) else unname(coef_means[nm])
      } else {
        cand <- c(paste0(col, ":", nm), paste0(nm, ":", col))
        hit  <- cand[cand %in% names(fe)]
        if (length(hit)) unname(fe[hit[1L]]) else 0
      }
    }, numeric(1L)),
    colnames(Xj)
  )

  cat(sprintf("  %-*s ~ %s\n", w, nm, hyper_rhs))
  print(gamma)
  cat("\n")
}

## ---------------------------------------------------------------------------
## 5. Empirical SD/variance of per-school coefficients vs lmer VarCorr
## ---------------------------------------------------------------------------
cat("--- Between-school SD of random coefficients vs lmer VarCorr ---\n")
vc <- as.data.frame(lme4::VarCorr(design$lmer_fit))
cat(sprintf("  %-16s  empirical_sd=%7.4f  empirical_var=%8.4f  lmer_sd=%7.4f  lmer_var=%8.4f\n",
            "(Intercept)",
            coef_sds["(Intercept)"], coef_vars["(Intercept)"],
            vc$sdcor[vc$var1 == "(Intercept)" & is.na(vc$var2)][1L],
            vc$vcov[vc$var1  == "(Intercept)" & is.na(vc$var2)][1L]))
for (nm in setdiff(design$re_coef_names, "(Intercept)")) {
  if (!nm %in% colnames(coef_df)) next
  lmer_row <- vc[vc$var1 == nm & is.na(vc$var2), ]
  cat(sprintf("  %-16s  empirical_sd=%7.4f  empirical_var=%8.4f  lmer_sd=%7.4f  lmer_var=%8.4f\n",
              nm,
              coef_sds[nm], coef_vars[nm],
              if (nrow(lmer_row)) lmer_row$sdcor[1L] else NA_real_,
              if (nrow(lmer_row)) lmer_row$vcov[1L]  else NA_real_))
}

## ---------------------------------------------------------------------------
## 6. lmer refitted on full-rank schools only (same subset as Prior_Setup)
## ---------------------------------------------------------------------------
cat("\n--- lmer refit: full-rank schools only ---\n")
full_rank_schools <- names(design$re_rank)[design$re_rank]
cat(sprintf("  Using %d of %d schools (dropping rank-deficient: %s)\n\n",
            length(full_rank_schools),
            nlevels(design$groups),
            paste(names(design$re_rank)[!design$re_rank], collapse = ", ")))

dat_fr <- subset(dat, school_id %in% full_rank_schools)
dat_fr$school_id <- droplevels(dat_fr$school_id)
fit_fr <- lme4::lmer(form_lmer, data = dat_fr, control = ctrl_bobyqa)
print(summary(fit_fr))

cat("\n--- VarCorr comparison: all schools vs full-rank schools only ---\n")
vc_fr <- as.data.frame(lme4::VarCorr(fit_fr))
cat(sprintf("  %-16s  all_schools_sd=%7.4f  full_rank_sd=%7.4f\n",
            "(Intercept)",
            vc$sdcor[vc$var1 == "(Intercept)" & is.na(vc$var2)][1L],
            vc_fr$sdcor[vc_fr$var1 == "(Intercept)" & is.na(vc_fr$var2)][1L]))
for (nm in setdiff(design$re_coef_names, "(Intercept)")) {
  row_all <- vc[vc$var1 == nm & is.na(vc$var2), ]
  row_fr  <- vc_fr[vc_fr$var1 == nm & is.na(vc_fr$var2), ]
  cat(sprintf("  %-16s  all_schools_sd=%7.4f  full_rank_sd=%7.4f\n",
              nm,
              if (nrow(row_all)) row_all$sdcor[1L] else NA_real_,
              if (nrow(row_fr))  row_fr$sdcor[1L]  else NA_real_))
}

## ===========================================================================
## 7. Optional stress test (NOT the lmerb example model): esl_observed RE
## ===========================================================================
cat("\n\n=== Section 7 (optional): esl_observed added as random slope ===\n\n")

dat_esl <- subset(
  dat,
  complete.cases(dat[, c("esl_observed")])
)
form_esl <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_a1 + distracted_ppvt +
  free_reduced_lunch:distracted_a1 + esl_observed +
  (1 + distracted_ppvt + distracted_a1 + esl_observed || school_id)

design_esl <- model_setup(form_esl, data = dat_esl, control = ctrl_bobyqa)
print(design_esl)

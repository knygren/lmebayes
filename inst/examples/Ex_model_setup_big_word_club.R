## model_setup() on bayesrules::big_word_club (private-school moderation)
##
## Two-level structure:
##   Level 1 (students):  y ~ b0[j] + b_age_c[j]*age_c
##   Level 2 (schools):   b0[j]      ~ private_school + title1 + free_reduced_lunch
##                        b_age_c[j] ~ 1

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Example requires the 'bayesrules' package.", call. = FALSE)
}

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat$age_c <- dat$age_months - mean(dat$age_months, na.rm = TRUE)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "age_c",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  (1 + age_c || school_id)

## ---------------------------------------------------------------------------
## 1. lmer fit: raw output
## ---------------------------------------------------------------------------
cat("--- lmer fit ---\n")
fit <- lme4::lmer(form_lmer, data = dat)
print(summary(fit))

cat("\n--- fixef(fit): population-level (gamma) estimates ---\n")
print(lme4::fixef(fit))

cat("\n--- coef(fit): per-group coefficients (fixef + ranef) ---\n")
print(coef(fit))

## ---------------------------------------------------------------------------
## 2. model_setup: structured view of the same model
## ---------------------------------------------------------------------------
design <- model_setup(form_lmer, data = dat)
print(design)

## ---------------------------------------------------------------------------
## 3. Random effects b[j]: first 10 schools
##    Columns match the RE predictors in the Measurement Model above
## ---------------------------------------------------------------------------
cat("--- Random effects b[j]: first 10", design$group_name, "---\n")
re_df <- as.data.frame(lme4::ranef(design$lmer_fit)[[design$group_name]])
print(utils::head(re_df, 10))

## ---------------------------------------------------------------------------
## 4. Gamma estimates organised to match the Random Effects Model above
##
##    Mapping:
##      intercept RE -- X_hyper column names map directly to fixef() names
##      slope RE     -- (Intercept) col -> fe[nm]  (population mean slope gamma_10)
##                      other cols      -> col:nm or nm:col interaction in fixef
## ---------------------------------------------------------------------------
cat("\n--- Random effects model (gamma estimates) ---\n")

fe          <- lme4::fixef(design$lmer_fit)
coef_df     <- coef(design$lmer_fit)[[design$group_name]]
coef_means  <- colMeans(coef_df)
coef_vars   <- apply(coef_df, 2L, var)
coef_sds    <- sqrt(coef_vars)
w           <- max(nchar(design$re_coef_names))

for (nm in design$re_coef_names) {
  Xj    <- design$X_hyper[[nm]]
  other <- setdiff(colnames(Xj), "(Intercept)")
  hyper_rhs <- if (length(other) == 0L) "1" else paste(c("1", other), collapse = " + ")

  gamma <- setNames(
    vapply(colnames(Xj), function(col) {
      if (nm == "(Intercept)") {
        if (col %in% names(fe)) unname(fe[col]) else 0
      } else if (col == "(Intercept)") {
        if (nm %in% names(fe)) unname(fe[nm])
        else if (nm %in% names(coef_means)) unname(coef_means[nm])
        else 0
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
##    For the intercept RE this is the between-school SD in mean scores.
##    For slope REs this is the between-school SD in the age/female effects.
## ---------------------------------------------------------------------------
cat("--- Between-school SD of random coefficients vs lmer VarCorr ---\n")
vc <- as.data.frame(lme4::VarCorr(design$lmer_fit))
cat(sprintf("  %-13s  empirical_sd=%7.4f  empirical_var=%8.4f  lmer_sd=%7.4f  lmer_var=%8.4f\n",
            "(Intercept)",
            coef_sds["(Intercept)"], coef_vars["(Intercept)"],
            vc$sdcor[vc$var1 == "(Intercept)" & is.na(vc$var2)][1L],
            vc$vcov[vc$var1  == "(Intercept)" & is.na(vc$var2)][1L]))
for (nm in setdiff(design$re_coef_names, "(Intercept)")) {
  if (!nm %in% colnames(coef_df)) next
  lmer_row <- vc[vc$var1 == nm & is.na(vc$var2), ]
  cat(sprintf("  %-13s  empirical_sd=%7.4f  empirical_var=%8.4f  lmer_sd=%7.4f  lmer_var=%8.4f\n",
              nm,
              coef_sds[nm], coef_vars[nm],
              if (nrow(lmer_row)) lmer_row$sdcor[1L] else NA_real_,
              if (nrow(lmer_row)) lmer_row$vcov[1L]  else NA_real_))
}

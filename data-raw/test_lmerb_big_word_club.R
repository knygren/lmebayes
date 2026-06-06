# model_setup() on bayesrules::big_word_club — three hyper matrices
#
# Design formula includes female:private_school (hyper / fixef calibration).
# vcov_re uses lmerb_default_vcov_formula (no fixed interaction; keeps ||).
#
#   score_ppvt ~ private_school + title1 + female:private_school +
#     (1 + age_c + female || school_id)
#
# Three X_hyper matrices (one row per school):
#   (Intercept) ~ 1 + private_school + title1  (47 x 3)
#   age_c         ~ 1                            (47 x 1)
#   female        ~ 1 + private_school           (47 x 2)

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules: install.packages('bayesrules')")
}
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload: install.packages('pkgload')")
}

pkgload::load_all(export_all = FALSE)

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$age_c <- dat$age_months - mean(dat$age_months, na.rm = TRUE)
dat$school_id <- factor(dat$school_id)

dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "age_c", "female",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)

n_schools <- nlevels(dat$school_id)
cat("big_word_club: n =", nrow(dat), ", schools =", n_schools, "\n")

form_lmer_interaction <- score_ppvt ~
  private_school + title1 +
  female:private_school +
  (1 + age_c + female || school_id)

# ---------------------------------------------------------------------------
# model_setup: three hyper design matrices
# ---------------------------------------------------------------------------
cat("\n--- model_setup: three X_hyper matrices ---\n")
design <- model_setup(form_lmer_interaction, data = dat)
print(design)

X_int <- design$X_hyper[["(Intercept)"]]
X_age <- design$X_hyper[["age_c"]]
X_fem <- design$X_hyper[["female"]]

stopifnot(
  inherits(design, "model_setup"),
  identical(design$re_coef_names, c("(Intercept)", "age_c", "female")),
  nrow(design$Z) == nrow(dat),
  ncol(design$Z) == 3L,
  identical(colnames(design$Z), design$re_coef_names),
  length(design$groups) == nrow(design$Z),
  length(design$y) == nrow(design$Z),
  isTRUE(all.equal(
    design$y,
    stats::model.response(lme4::lFormula(form_lmer_interaction, data = dat)$fr)
  )),
  design$Z[, "(Intercept)"] == 1,
  isTRUE(all.equal(design$Z[, "age_c"], dat$age_c)),
  isTRUE(all.equal(design$Z[, "female"], dat$female)),
  nrow(X_int) == n_schools,
  ncol(X_int) == 3L,
  all(c("(Intercept)", "private_school", "title1") %in% colnames(X_int)),
  nrow(X_age) == n_schools,
  ncol(X_age) == 1L,
  colnames(X_age) == "(Intercept)",
  nrow(X_fem) == n_schools,
  ncol(X_fem) == 2L,
  all(c("(Intercept)", "private_school") %in% colnames(X_fem))
)

cat("\nIntercept hyper (head):\n")
print(utils::head(X_int))
cat("\nage_c hyper (head):\n")
print(utils::head(X_age))
cat("\nfemale hyper (head):\n")
print(utils::head(X_fem))

schools_5 <- levels(design$groups)[1:5]
idx5 <- design$groups %in% schools_5
cat("\n========== y and Z: first 5 schools (",
    paste(schools_5, collapse = ", "), ", n = ", sum(idx5), ") ==========\n", sep = "")
print(cbind(y = design$y[idx5], design$Z[idx5, , drop = FALSE]))

sid2 <- levels(design$groups)[2L]
idx2 <- design$groups == sid2
cat("\n========== y and Z: second school only (school_id =", sid2,
    ", n =", sum(idx2), ") ==========\n")
print(cbind(y = design$y[idx2], design$Z[idx2, , drop = FALSE]))

stopifnot(
  inherits(design$lmer_fit, "lmerMod"),
  !is.null(design$varcorr),
  length(design$vcov_re) == 3L,
  identical(names(design$vcov_re), design$re_coef_names),
  is.numeric(design$residual_var),
  design$residual_var > 0,
  grepl("||", deparse1(design$vcov_formula), fixed = TRUE)
)

cat("\nvcov_formula (RE calibration, same || as design):\n")
print(design$vcov_formula)

cat("\n--- lmer variance components (from model_setup) ---\n")
print(design$varcorr)

# ---------------------------------------------------------------------------
# lmer reference: gamma vectors align with X_hyper columns
# ---------------------------------------------------------------------------
cat("\n--- lmer: hyper coefficients ---\n")
fe <- lme4::fixef(design$lmer_fit)
cat("fixef:\n")
print(fe)

gamma_int <- fe[colnames(X_int)]
# lmer has private_school:female only (no female intercept in fixef) -> gamma_0 = 0
gamma_fem <- c("(Intercept)" = 0, private_school = unname(fe["private_school:female"]))
mu_int_school1 <- as.numeric(X_int[1L, , drop = TRUE] %*% gamma_int)
mu_fem_school1 <- as.numeric(X_fem[1L, , drop = TRUE] %*% gamma_fem)

cat("\nSchool 1 prior means from X_hyper %*% gamma:\n")
cat("  mu_(Intercept) =", mu_int_school1, "\n")
cat("  mu_age_c       = 0  (intercept-only X, no gamma in lmer fixef)\n")
cat("  mu_female      =", mu_fem_school1,
    "  (gamma_fem: Intercept=0, private_school=", gamma_fem["private_school"], ")\n", sep = "")

cat("\nmodel_setup (big_word_club): OK\n")

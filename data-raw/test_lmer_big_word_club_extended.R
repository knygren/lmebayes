# Extended big_word_club lmer — one female:Z interaction; age/female RE only.
#
# Block likelihood: y ~ 1 + age_c + female
# Hyper: intercept ~ Z (X_nbhd); female ~ private_school in prior mean (not in likelihood X)
#
#   score_ppvt ~
#     private_school + title1 + free_reduced_lunch + female:private_school +
#     (1 + age_c + female || school_id)
#
# No fixed age_c or female mains — student slopes enter only via random effects,
# except female moderation by private_school (the single X:Z interaction).

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

cat("big_word_club extended model: n =", nrow(dat),
    ", schools =", nlevels(dat$school_id), "\n")

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  female:private_school +
  (1 + age_c + female || school_id)

stopifnot(lmebayes:::is_single_factor_model(form_lmer, data = dat))

lme4_comps <- lmebayes:::get_lme4_components(form_lmer, data = dat)


Print(lme4_comps$Z_random_sparse)

X_nbhd <- lmebayes:::extract_lme4_fixed_group_matrix(lme4_comps, "school_id")
cat("\n--- X_nbhd (school-constant fixed columns only) ---\n")
cat("dim:", paste(dim(X_nbhd), collapse = " x "), "\n")
print(colnames(X_nbhd))

# Level-1 rank: student likelihood X only (private_school is hyper, not in X).
form_block <- score_ppvt ~ age_c + female
mf_block <- stats::model.frame(form_block, data = dat)
x_likelihood <- stats::model.matrix(attr(mf_block, "terms"), mf_block)
grp <- lme4_comps$groups$school_id
stopifnot(nrow(x_likelihood) == length(grp))
cat("\nLikelihood design (block_lmb X), cols:",
    paste(colnames(x_likelihood), collapse = ", "), "\n")

id_check <- block_check_identifiability_xy(
  x = x_likelihood,
  block = grp,
  X_nbhd = X_nbhd,
  on_failure = "warn"
)

fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)

if (lme4::isSingular(fit_lmer)) {
  message("lmer fit is singular — check VarCorr; RE variances may be on boundary.")
}

cat("\n--- lmer: extended big_word_club ---\n")
print(summary(fit_lmer))

cat("\nFixed effects:\n")
print(lme4::fixef(fit_lmer))

cat("\nRandom effects (head):\n")
print(utils::head(lme4::ranef(fit_lmer)$school_id))

cat("\nVariance components:\n")
print(lme4::VarCorr(fit_lmer))

example_school <- rownames(X_nbhd)[1L]
slice <- lmebayes:::extract_lme4_submatrices(
  lme4_comps, "school_id", target_level = example_school
)
cat("\n--- Example school", example_school, "slice ---\n")
cat("X_fixed_level1 cols:", paste(colnames(slice$X_fixed_level1), collapse = ", "), "\n")
cat("X_fixed_level2 cols:", paste(colnames(slice$X_fixed_level2), collapse = ", "), "\n")
cat("Z_random_level1 dim:", paste(dim(slice$Z_random_level1), collapse = " x "), "\n")

cat("\nlmer extended (big_word_club): OK\n")

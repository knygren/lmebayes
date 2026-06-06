# Reference lme4 fit for bayesrules::big_word_club -- intercept-as-outcomes
# with one female:Z interaction and glmerb_utilities extraction.
#
# Block / Bayesian layers (two formulas):
#   Likelihood (per school):  y ~ 1 + age_c + female
#   Hyper:  intercept ~ private_school + title1 + free_reduced_lunch  (X_nbhd)
#           age_c ~ 1
#           female ~ private_school  (prior mean of female coef, NOT in likelihood X)
#
# lmer reference fit (single combined formula -- interaction is a fixed column):
#   score_ppvt ~
#     private_school + title1 + free_reduced_lunch + female:private_school +
#     (1 + age_c + female || school_id)
#
# Uses glmerb_utilities to extract X_nbhd (1 ~ Z), per-school student X,
# school attributes, and random-effects Z.

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

cat("big_word_club analysis n =", nrow(dat),
    ", schools =", nlevels(dat$school_id), "\n")

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  female:private_school +
  (1 + age_c + female || school_id)

stopifnot(lmebayes:::is_single_factor_model(form_lmer, data = dat))

lme4_comps <- get_lme4_components(form_lmer, data = dat)

cat("\n--- Z_random labels: View(lme4_comps$Z_random_column_map) ---\n")
cat("              also: View(lme4_comps$Z_random_row_map)\n")
if (interactive()) {
  View(lme4_comps$Z_random_column_map)
  View(lme4_comps$Z_random_row_map)
} else {
  print(head(lme4_comps$Z_random_column_map))
}

X_nbhd <- lmebayes:::extract_lme4_fixed_group_matrix(lme4_comps, "school_id")
cat("\n--- X_nbhd (1 ~ Z), dim =", paste(dim(X_nbhd), collapse = " x "), "---\n")
print(head(X_nbhd))

example_school <- rownames(X_nbhd)[1L]
slice <- lmebayes:::extract_lme4_submatrices(
  lme4_comps, "school_id", target_level = example_school
)
cat("\n--- Example school", example_school, "---\n")
cat("n_observations:", slice$n_observations, "\n")
cat("X_fixed_level1 cols:", paste(colnames(slice$X_fixed_level1), collapse = ", "), "\n")
cat("X_fixed_level2:\n")
print(slice$X_fixed_level2)
cat("Z_random_level1 dim:", paste(dim(slice$Z_random_level1), collapse = " x "), "\n")

# Likelihood X: student terms only. private_school moderates female in the prior
# (mu_b["female"] = private_school_j * gamma_ps), not as female:private_school.
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
  message("lmer fit is singular -- check VarCorr; RE variances may be on boundary.")
}

cat("\n--- lmer: big_word_club intercept-as-outcomes ---\n")
print(summary(fit_lmer))

fe <- lme4::fixef(fit_lmer)
re_school <- lme4::ranef(fit_lmer)$school_id

cat("\nFixed effects:\n")
print(fe)

cat("\nRandom effects (per-school deviations, head):\n")
print(utils::head(re_school))

vc_df <- as.data.frame(lme4::VarCorr(fit_lmer))
var_re_int <- vc_df$vcov[
  vc_df$grp == "school_id" & vc_df$var1 == "(Intercept)"
]
var_re_age <- vc_df$vcov[
  vc_df$grp == "school_id.1" & vc_df$var1 == "age_c"
]
if (length(var_re_age) == 0L) {
  var_re_age <- vc_df$vcov[
    vc_df$grp == "school_id" & vc_df$var1 == "age_c"
  ]
}
var_re_female <- vc_df$vcov[
  vc_df$grp == "school_id.2" & vc_df$var1 == "female"
]
if (length(var_re_female) == 0L) {
  var_re_female <- vc_df$vcov[
    vc_df$grp == "school_id" & vc_df$var1 == "female"
  ]
}
dispersion_lmer <- vc_df$vcov[vc_df$grp == "Residual"]

if (var_re_age <= 0) {
  message("Zero variance for age_c RE; using small prior variance in block_lmb.")
  var_re_age <- max(var_re_age, 1e-4)
}
if (length(var_re_female) == 0L || var_re_female <= 0) {
  message("Zero variance for female RE; using small prior variance in block_lmb.")
  var_re_female <- max(var_re_female, 1e-6)
}

z_coef_names <- c(
  "(Intercept)", "private_school", "title1", "free_reduced_lunch"
)
stopifnot(all(z_coef_names %in% names(fe)))

int_name <- grep("private_school.*female|female.*private_school", names(fe), value = TRUE)
stopifnot(length(int_name) == 1L)
gamma_ps <- fe[int_name]

gamma_z <- fe[z_coef_names]
mu_by_school <- X_nbhd %*% gamma_z

coef_names_block <- colnames(x_likelihood)
stopifnot(all(c("(Intercept)", "age_c", "female") %in% coef_names_block))

pfamily_list <- lapply(rownames(X_nbhd), function(sid) {
  ps <- X_nbhd[sid, "private_school"]
  mu_b <- setNames(
    c(mu_by_school[sid, 1L], 0, ps * gamma_ps),
    coef_names_block
  )
  Sigma_b <- setNames(
    c(var_re_int, var_re_age, var_re_female),
    coef_names_block
  )
  Sigma_b <- diag(as.numeric(Sigma_b))
  dimnames(Sigma_b) <- list(coef_names_block, coef_names_block)
  glmbayes::dNormal(mu = mu_b, Sigma = Sigma_b, dispersion = dispersion_lmer)
})
names(pfamily_list) <- rownames(X_nbhd)

set.seed(42)
n_draw <- 1000L

out_blmb <- block_lmb(
  form_block,
  block = "school_id",
  pfamily_list = pfamily_list,
  data = dat,
  n = n_draw,
  use_parallel = FALSE
)

stopifnot(inherits(out_blmb, "blmb"), length(out_blmb) == nrow(X_nbhd))

cm <- lmebayes:::.blmb_coef_means_matrix(out_blmb)
ps_by_school <- X_nbhd[rownames(re_school), "private_school"]
coef_lmer_school <- cbind(
  as.numeric(mu_by_school) + re_school[, "(Intercept)"],
  re_school[, "age_c"],
  ps_by_school * gamma_ps + re_school[, "female"]
)
colnames(coef_lmer_school) <- coef_names_block
rownames(coef_lmer_school) <- rownames(re_school)

cat("\n--- block_lmb posterior coefficient means (head) ---\n")
print(utils::head(cm))

cat("\n--- lmer coefficients by school (intercept/age BLUP; female = ps*gamma_ps + RE, head) ---\n")
print(utils::head(coef_lmer_school))

diff_mat <- cm[rownames(coef_lmer_school), coef_names_block, drop = FALSE] -
  coef_lmer_school
cat("\n--- max |block_lmb mean - lmer BLUP| ---\n")
print(max(abs(diff_mat)))

cat("\nlmer + glmerb_utilities (big_word_club): OK\n")

# Reference lme4 fit for sleepstudy — uncorrelated random intercept and random
# Days slope by Subject (diagonal RE, matching test_lmer_iris.R).
#
# Days is centered at the sample mean before fitting so posterior subject
# coefficients are less correlated (intercept = reaction at mean day).
#
# Canonical lme4 example: Reaction ~ Days + (Days | Subject)
# This script uses uncorrelated RE so block_lmb priors use diagonal Sigma:
#   Reaction ~ Days_c + (1 + Days_c || Subject)
#
# block_lmb: Reaction ~ Days_c, BY Subject, dNormal priors from lmer VarCorr

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}

data(sleepstudy, package = "lme4")

days_mean <- mean(sleepstudy$Days)
sleepstudy$Days_c <- sleepstudy$Days - days_mean
cat(
  "Centered Days at mean =", days_mean,
  "; Days_c range:", paste(range(sleepstudy$Days_c), collapse = " to "), "\n"
)

form_lmer <- Reaction ~ Days_c + (1 + Days_c || Subject)

fit_lmer <- lme4::lmer(
  form_lmer,
  data = sleepstudy,
  REML = TRUE
)

cat("\n--- lmer: Days_c + (1 + Days_c || Subject) ---\n")
print(summary(fit_lmer))

cat("\nFixed effects:\n")
print(lme4::fixef(fit_lmer))

cat("\nRandom effects (per-Subject deviations):\n")
print(lme4::ranef(fit_lmer)$Subject)

cat("\nVariance components:\n")
print(lme4::VarCorr(fit_lmer))

# ---------------------------------------------------------------------------
# Smoke checks
# ---------------------------------------------------------------------------
stopifnot(inherits(fit_lmer, "lmerMod"))
stopifnot(!lme4::isSingular(fit_lmer))
stopifnot(nrow(sleepstudy) == 180L)

fe <- lme4::fixef(fit_lmer)
stopifnot(length(fe) == 2L)
stopifnot(all(c("(Intercept)", "Days_c") %in% names(fe)))

re <- lme4::ranef(fit_lmer)
stopifnot(length(re) == 1L)
stopifnot("Subject" %in% names(re))
re_sub <- re$Subject
stopifnot(nrow(re_sub) == 18L, ncol(re_sub) == 2L)
stopifnot(all(c("(Intercept)", "Days_c") %in% colnames(re_sub)))

vc <- lme4::VarCorr(fit_lmer)
stopifnot(length(vc) == 2L) # Subject + Subject.1

cat("\nlmer (sleepstudy, centered Days, uncorrelated RE): OK\n")

# ---------------------------------------------------------------------------
# block_lmb with dNormal priors calibrated from lmer (BY Subject)
# ---------------------------------------------------------------------------
# Per-block prior:
#   mu         = lmer fixed effects (population intercept + Days_c slope)
#   Sigma      = diag(RE variance components for intercept and Days_c slope)
#   dispersion = (residual Std.Dev.)^2

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload: install.packages('pkgload')")
}
pkgload::load_all(export_all = FALSE)

vc_df <- as.data.frame(lme4::VarCorr(fit_lmer))
var_re_int <- vc_df$vcov[
  vc_df$grp == "Subject" & vc_df$var1 == "(Intercept)"
]
var_re_slope <- vc_df$vcov[
  vc_df$grp == "Subject.1" & vc_df$var1 == "Days_c"
]
dispersion_lmer <- vc_df$vcov[vc_df$grp == "Residual"]

coef_names <- names(fe)
mu_lmer <- as.numeric(fe)
names(mu_lmer) <- coef_names

Sigma_lmer <- diag(c(var_re_int, var_re_slope))
dimnames(Sigma_lmer) <- list(coef_names, coef_names)

cat("\n--- block_lmb priors from lmer ---\n")
cat("mu:\n")
print(mu_lmer)
cat("\nSigma (diagonal, RE variances):\n")
print(Sigma_lmer)
cat("\ndispersion (residual variance):", dispersion_lmer, "\n")

form_block <- Reaction ~ Days_c

pfamily_list <- lapply(rownames(re_sub), function(.x) {
  glmbayes::dNormal(
    mu = mu_lmer,
    Sigma = Sigma_lmer,
    dispersion = dispersion_lmer
  )
})
names(pfamily_list) <- rownames(re_sub)

set.seed(42)
n_draw <- 1000L

out_blmb <- lmbBlock(
  form_block,
  block = "Subject",
  pfamily_list = pfamily_list,
  data = sleepstudy,
  n = n_draw,
  use_parallel = FALSE
)

stopifnot(inherits(out_blmb, "blmb"), length(out_blmb) == 18L)
stopifnot(inherits(out_blmb[[1L]], "lmb"))
stopifnot(nrow(out_blmb[[1L]]$coefficients) == n_draw)

cm <- attr(summary(out_blmb), "coef_means")
stopifnot(nrow(cm) == 18L, ncol(cm) == 2L)
stopifnot(all(colnames(cm) == coef_names))

cat("\n--- block_lmb posterior coefficient means (by Subject) ---\n")
print(cm)

coef_lmer_by_subject <- re_sub + matrix(
  fe,
  nrow = nrow(re_sub),
  ncol = length(fe),
  byrow = TRUE
)

cat("\n--- lmer coefficients by Subject (fixef + ranef) ---\n")
print(coef_lmer_by_subject)

diff_mat <- cm[rownames(coef_lmer_by_subject), , drop = FALSE] - coef_lmer_by_subject
cat("\n--- max |block_lmb mean - lmer BLUP| ---\n")
print(max(abs(diff_mat)))

cat("\nblock_lmb (sleepstudy, centered Days, dNormal priors from lmer): OK\n")

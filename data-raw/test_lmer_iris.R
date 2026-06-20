# Reference lme4 fit for iris — uncorrelated random intercept and random
# Sepal.Width slope by Species (Petal.Length omitted).
#
# Contrast with data-raw/test_block_lmb_iris.R:
#   block_lmb: Sepal.Length ~ Sepal.Width + Petal.Length, BY Species
#              (3 separate regressions, no cross-species pooling)
#   lmer:      Sepal.Length ~ Sepal.Width + (1 | Species) +
#                (0 + Sepal.Width | Species)
#              (species-specific intercept and Sepal.Width slope; RE
#               variances uncorrelated — identifiable with g = 3 groups)
#
# Equivalent shorthand: Sepal.Length ~ Sepal.Width + (1 + Sepal.Width || Species)

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}

data("iris", package = "datasets")

form_lmer <- Sepal.Length ~ Sepal.Width + (1 + Sepal.Width || Species)

fit_lmer <- lme4::lmer(
  form_lmer,
  data = iris,
  REML = TRUE
)

cat("\n--- lmer: Sepal.Width + (1 + Sepal.Width || Species) ---\n")
print(summary(fit_lmer))

cat("\nFixed effects:\n")
print(lme4::fixef(fit_lmer))

cat("\nRandom effects (per-species deviations):\n")
print(lme4::ranef(fit_lmer))

cat("\nVariance components:\n")
print(lme4::VarCorr(fit_lmer))

# ---------------------------------------------------------------------------
# Smoke checks (parallel structure to test_block_lmb_iris.R)
# ---------------------------------------------------------------------------
stopifnot(inherits(fit_lmer, "lmerMod"))
stopifnot(!lme4::isSingular(fit_lmer))
stopifnot(nrow(iris) == 150L)

fe <- lme4::fixef(fit_lmer)
stopifnot(length(fe) == 2L)
stopifnot(all(c("(Intercept)", "Sepal.Width") %in% names(fe)))

re <- lme4::ranef(fit_lmer)
stopifnot(length(re) == 1L)
stopifnot("Species" %in% names(re))
re_sp <- re$Species
stopifnot(nrow(re_sp) == 3L, ncol(re_sp) == 2L)
stopifnot(all(c("(Intercept)", "Sepal.Width") %in% colnames(re_sp)))
stopifnot(all(rownames(re_sp) == levels(iris$Species)))

vc <- lme4::VarCorr(fit_lmer)
stopifnot(length(vc) == 2L) # Species (intercept) + Species.1 (slope)

cat("\nlmer (iris, uncorrelated random intercept + Sepal.Width slope): OK\n")

# ---------------------------------------------------------------------------
# block_lmb with dNormal priors calibrated from lmer (BY Species)
# ---------------------------------------------------------------------------
# Same response and fixed predictor as lmer; per-block prior:
#   mu         = lmer fixed effects (population intercept + Sepal.Width slope)
#   Sigma      = diag(RE variance components for intercept and Sepal.Width slope)
#   dispersion = (residual Std.Dev.)^2

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload: install.packages('pkgload')")
}
pkgload::load_all(export_all = FALSE)

vc_df <- as.data.frame(lme4::VarCorr(fit_lmer))
var_re_int <- vc_df$vcov[
  vc_df$grp == "Species" & vc_df$var1 == "(Intercept)"
]
var_re_slope <- vc_df$vcov[
  vc_df$grp == "Species.1" & vc_df$var1 == "Sepal.Width"
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

form_block <- Sepal.Length ~ Sepal.Width

pfamily_list <- lapply(levels(iris$Species), function(.x) {
  glmbayes::dNormal(
    mu = mu_lmer,
    Sigma = Sigma_lmer,
    dispersion = dispersion_lmer
  )
})
names(pfamily_list) <- levels(iris$Species)

set.seed(42)
n_draw <- 50L

out_blmb <- lmbBlock(
  form_block,
  block = "Species",
  pfamily_list = pfamily_list,
  data = iris,
  n = n_draw,
  use_parallel = FALSE
)

stopifnot(inherits(out_blmb, "blmb"), length(out_blmb) == 3L)
stopifnot(inherits(out_blmb[[1L]], "lmb"))
stopifnot(nrow(out_blmb[[1L]]$coefficients) == n_draw)

cm <- attr(summary(out_blmb), "coef_means")
stopifnot(nrow(cm) == 3L, ncol(cm) == 2L)
stopifnot(all(colnames(cm) == coef_names))

cat("\n--- block_lmb posterior coefficient means (by Species) ---\n")
print(cm)

# lmer species-specific coefficients: fixed effects + random effects (BLUPs)
coef_lmer_by_species <- re_sp + matrix(
  fe,
  nrow = nrow(re_sp),
  ncol = length(fe),
  byrow = TRUE
)

cat("\n--- lmer coefficients by Species (fixef + ranef) ---\n")
print(coef_lmer_by_species)

cat("\nblock_lmb (iris, dNormal priors from lmer): OK\n")

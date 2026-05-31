# =============================================================================
# Test: rindepNormalGamma_reg_multi on iris
# Model: y = cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)
#        x = ~(Intercept) + versicolor + virginica  (p = 3, l1 = 4)
# n = 1000 iid draws; compare posterior means to lm() MLE
# =============================================================================

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
pkgload::load_all(export_all = FALSE)

set.seed(42)
n_draw <- 1000L

# --- Data -------------------------------------------------------------------
data(iris)

# Response: all four numeric columns
y <- as.matrix(iris[, 1:4])
colnames(y) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")

# Design: intercept + Species dummies (setosa is baseline)  p = 3
x <- model.matrix(~ Species, data = iris)
cat("y:", nrow(y), "x", ncol(y), "(", paste(colnames(y), collapse = ", "), ")\n")
cat("x:", nrow(x), "x", ncol(x), "(", paste(colnames(x), collapse = ", "), ")\n\n")

# --- OLS reference ----------------------------------------------------------
lm_fit <- lm(y ~ Species, data = iris)
lm_coef <- coef(lm_fit)    # p x l1 matrix (predictors x responses)
cat("=== OLS (lm) coefficients ===\n")
print(round(t(lm_coef), 4))

lm_sigma2 <- sapply(summary(lm_fit), function(s) s$sigma^2)
cat("\n=== OLS residual variance (sigma^2) per response ===\n")
print(round(lm_sigma2, 4))

# --- Prior via Prior_Setup (one per response column) ------------------------
x_df <- iris["Species"]

prior_list <- lapply(seq_len(ncol(y)), function(j) {
  df_j <- cbind(x_df, y_j = y[, j])
  ps <- Prior_Setup(y_j ~ Species, family = gaussian(), data = df_j)
  list(
    mu    = as.numeric(ps$mu),
    Sigma = ps$Sigma,
    shape = ps$shape_ING,
    rate  = ps$rate
  )
})

cat("\n=== Prior means per response column ===\n")
prior_mu_mat <- do.call(rbind, lapply(prior_list, `[[`, "mu"))
rownames(prior_mu_mat) <- colnames(y)
colnames(prior_mu_mat) <- colnames(x)
print(round(prior_mu_mat, 4))

# --- Draw 1000 iid samples --------------------------------------------------
cat("\n=== Drawing", n_draw, "iid samples ===\n")
out <- rindepNormalGamma_reg_multi(
  n          = n_draw,
  y          = y,
  x          = x,
  prior_list = prior_list,
  family     = gaussian(),
  use_parallel = FALSE,
  progbar    = FALSE
)

# --- Check structure --------------------------------------------------------
stopifnot(inherits(out, "mrglmb"))
stopifnot(length(out) == ncol(y))
for (j in seq_len(ncol(y))) {
  fit <- out[[j]]
  stopifnot(inherits(fit, "rglmb"))
  M <- fit$coefficients
  stopifnot(is.matrix(M), nrow(M) == n_draw, ncol(M) == ncol(x))
}
cat("Structure checks passed.\n")

# --- Posterior means --------------------------------------------------------
post_mean <- do.call(rbind, lapply(out, function(fit) colMeans(fit$coefficients)))
rownames(post_mean) <- colnames(y)
colnames(post_mean) <- colnames(x)

print(round(t(lm_coef), 4))

cat("\n=== Posterior means (", n_draw, "draws) ===\n")
print(round(post_mean, 4))

cat("\n=== Difference: posterior mean - OLS ===\n")
diff_mat <- post_mean - t(lm_coef)
print(round(diff_mat, 4))

# --- Posterior mean dispersion ----------------------------------------------
disp_post_mean <- vapply(out, function(fit) mean(fit$dispersion), numeric(1))
names(disp_post_mean) <- colnames(y)
cat("\n=== Posterior mean dispersion (sigma^2) per response ===\n")
print(round(disp_post_mean, 4))

cat("\n=== OLS residual variance (sigma^2) per response ===\n")
print(round(lm_sigma2, 4))

cat("\nDone.\n")

# Smoke test: rNormal_reg_multi on iris (n = 200, l1 = 4, p = 3)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
pkgload::load_all(export_all = FALSE)

set.seed(42)
n_draw <- 200L

y <- as.matrix(iris[, 1:4])
colnames(y) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
x <- model.matrix(~ Species, data = iris)
x_df <- iris["Species"]

prior_list <- lapply(seq_len(ncol(y)), function(j) {
  df_j <- cbind(x_df, y_j = y[, j])
  ps   <- Prior_Setup(y_j ~ Species, family = gaussian(), data = df_j)
  list(mu = as.numeric(ps$mu), Sigma = ps$Sigma, dispersion = ps$dispersion)
})

out <- rNormal_reg_multi(
  n = n_draw, y = y, x = x, prior_list = prior_list,
  family = gaussian(), use_parallel = FALSE, progbar = FALSE
)

stopifnot(inherits(out, "mrglmb"), length(out) == 4L)
for (j in seq_len(ncol(y))) {
  stopifnot(inherits(out[[j]], "rglmb"))
  stopifnot(nrow(out[[j]]$coefficients) == n_draw, ncol(out[[j]]$coefficients) == 3L)
}
cat("rNormal_reg_multi: OK\n")

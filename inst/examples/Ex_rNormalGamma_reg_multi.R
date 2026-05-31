## rNormalGamma_reg_multi: multivariate response (iris)
## y: four numeric columns; x: intercept + Species (p = 3)

set.seed(42)

y <- as.matrix(iris[, 1:4])
colnames(y) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
x <- model.matrix(~ Species, data = iris)

x_df <- iris["Species"]
prior_list <- lapply(seq_len(ncol(y)), function(j) {
  df_j <- cbind(x_df, y_j = y[, j])
  ps   <- Prior_Setup(y_j ~ Species, family = gaussian(), data = df_j)
  list(
    mu = as.numeric(ps$mu),
    Sigma = ps$Sigma_0,
    shape = ps$shape,
    rate = ps$rate
  )
})

out <- rNormalGamma_reg_multi(
  n = 500, y = y, x = x, prior_list = prior_list,
  family = gaussian(), use_parallel = FALSE, progbar = FALSE
)

class(out)
names(out)
summary(out)

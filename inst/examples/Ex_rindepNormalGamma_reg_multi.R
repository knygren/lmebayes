## rindepNormalGamma_reg_multi: multivariate response example
## Data: iris  (n = 150)
## y: all four numeric measurements (l1 = 4 response columns)
## x: intercept + Species dummies (p = 3 predictors, shared across all columns)
## Returns an "mrglmb" object: a named list of rglmb fits, one per response column.

set.seed(42)

## --- Prepare data -----------------------------------------------------------
y <- as.matrix(iris[, 1:4])
colnames(y) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
x <- model.matrix(~ Species, data = iris)   # 150 x 3

## --- OLS reference ----------------------------------------------------------
lm_fit  <- lm(y ~ Species, data = iris)
lm_coef <- coef(lm_fit)   # p x l1 matrix

summary(lm_fit)

## --- One prior_list per response column (via Prior_Setup) -------------------
x_df <- iris["Species"]
prior_list <- lapply(seq_len(ncol(y)), function(j) {
  df_j <- cbind(x_df, y_j = y[, j])
  ps   <- Prior_Setup(y_j ~ Species, family = gaussian(), data = df_j)
  list(mu = as.numeric(ps$mu), Sigma = ps$Sigma,
       shape = ps$shape_ING, rate = ps$rate)
})

## --- Draw 1000 iid samples --------------------------------------------------
out <- rindepNormalGamma_reg_multi(
  n = 1000, y = y, x = x, prior_list = prior_list,
  family = gaussian(), use_parallel = FALSE, progbar = FALSE
)

## out is an "mrglmb" object: a named list of rglmb fits
class(out)          # "mrglmb"
names(out)          # "Sepal.Length" "Sepal.Width" "Petal.Length" "Petal.Width"

## Access a single response fit (class "rglmb")
sl_fit <- out[["Sepal.Length"]]
dim(sl_fit$coefficients)   # 1000 x 3  (draws x predictors)
length(sl_fit$dispersion)  # 1000

## --- Posterior means of coefficients (all responses) -----------------------
post_mean <- do.call(rbind, lapply(out, function(fit) colMeans(fit$coefficients)))
colnames(post_mean) <- colnames(x)
round(post_mean, 3)
round(t(lm_coef), 3)   # OLS reference

## --- Posterior mean dispersion per response ---------------------------------
disp_post <- vapply(out, function(fit) mean(fit$dispersion), numeric(1))
round(disp_post, 4)

summary(out)

## summary.rGamma_reg: dGamma prior (dispersion-only; coefficients fixed)
## rGamma_reg, rglmb, and rlmb use summary.rGamma_reg when prior is dGamma.
data("Boston", package = "MASS")

predictors <- setdiff(names(Boston), "medv")
Boston_centered <- Boston
Boston_centered[predictors] <- scale(Boston[predictors], center = TRUE, scale = FALSE)

form <- medv ~
  crim + zn +
  indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

ps.boston <- glmbayes::Prior_Setup(form, gaussian(), data = Boston_centered)
rate_dg <- if (!is.null(ps.boston$rate_gamma)) ps.boston$rate_gamma else ps.boston$rate

y <- ps.boston$y
x <- as.matrix(ps.boston$x)
wt <- rep(1, length(y))

out1 <- rGamma_reg(
  n = 200,
  y = y,
  x = x,
  prior_list = list(beta = ps.boston$coefficients, shape = ps.boston$shape, rate = rate_dg),
  offset = rep(0, length(y)),
  weights = wt,
  family = gaussian()
)
summary(out1)

out2 <- glmbayes::rglmb(
  n = 200,
  y = y,
  x = x,
  pfamily = glmbayes::dGamma(shape = ps.boston$shape, rate = rate_dg, beta = ps.boston$coefficients),
  weights = wt,
  family = gaussian()
)
summary(out2)

out3 <- glmbayes::rlmb(
  n = 200,
  y = y,
  x = x,
  pfamily = glmbayes::dGamma(shape = ps.boston$shape, rate = rate_dg, beta = ps.boston$coefficients),
  weights = wt
)
summary(out3)

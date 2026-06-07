## Row-block lmb (iris): BY Species — SAS-style separate regressions

set.seed(42)
data("iris", package = "datasets")

ps_block <- block_prior_setup(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  data = iris,
  family = gaussian()
)

pfamily_list <- lapply(ps_block, function(ps) {
  glmbayesCore::dNormal_Gamma(
    mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate
  )
})

out_blmb <- block_lmb(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  pfamily_list = pfamily_list,
  data = iris,
  n = 150L,
  use_parallel = FALSE
)

print(out_blmb)
summary(out_blmb)

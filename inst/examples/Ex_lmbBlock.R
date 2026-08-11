## Row-block lmb (iris): BY Species — SAS-style separate regressions

set.seed(42)
data("iris", package = "datasets")

ps_block <- Prior_SetupGroup(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  group = "Species",
  data = iris,
  family = gaussian()
)

out_blmb <- lmbBlock(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  pfamily_list = pfamily_list(ps_block),
  data = iris,
  n = 150L,
  use_parallel = FALSE
)

print(out_blmb)
summary(out_blmb)

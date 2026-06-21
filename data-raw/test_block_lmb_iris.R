# Smoke test: lmbBlock on iris (BY Species)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
pkgload::load_all(export_all = FALSE)

data("iris", package = "datasets")

set.seed(42)
n_draw <- 50L

ps_block <- Prior_SetupBlock(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  data = iris,
  family = gaussian()
)
stopifnot(inherits(ps_block, "block_PriorSetup"), length(ps_block) == 3L)

pfamily_list <- lapply(ps_block, function(ps) {
  dNormal_Gamma(
    mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate
  )
})

out <- lmbBlock(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  block = "Species",
  pfamily_list = pfamily_list,
  data = iris,
  n = n_draw,
  use_parallel = FALSE
)

stopifnot(inherits(out, "blmb"), length(out) == 3L)
stopifnot(inherits(out[[1L]], "lmb"))
stopifnot(nrow(out[[1L]]$coefficients) == n_draw)

cm <- attr(summary(out), "coef_means")
stopifnot(nrow(cm) == 3L, ncol(cm) == 3L)

print(out)

# block_rNormalGLM is glmbayesCore-only (not re-exported from lmebayes)
stopifnot(exists("block_rNormalGLM", where = asNamespace("glmbayesCore")))

cat("lmbBlock (iris): OK\n")

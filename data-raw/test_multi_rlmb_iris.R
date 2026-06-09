# Smoke tests: multi_* functions on iris (n = 100, l1 = 4, p = 3)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
if (!requireNamespace("glmbayesCore", quietly = TRUE)) {
  stop("Install glmbayesCore.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

set.seed(42)
n_draw <- 100L

ps_multi <- glmbayesCore::multi_prior_setup(
  cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species,
  data = iris,
  family = gaussian()
)
stopifnot(inherits(ps_multi, "multi_PriorSetup"), length(ps_multi) == 4L)

x <- ps_multi[[1L]]$x
y <- do.call(cbind, lapply(ps_multi, function(ps) ps$y))
colnames(y) <- names(ps_multi)

common <- list(
  n = n_draw, y = y, x = x,
  family = gaussian(), use_parallel = FALSE, progbar = FALSE
)

pfamily_list <- lapply(ps_multi, function(ps) {
  glmbayesCore::dNormal_Gamma(
    mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate
  )
})
out_rlmb <- multi_rlmb(
  n = n_draw, y = y, x = x, pfamily_list = pfamily_list,
  use_parallel = FALSE, progbar = FALSE
)
stopifnot(inherits(out_rlmb, "mrglmb"), length(out_rlmb) == 4L)
stopifnot(inherits(out_rlmb[[1L]], "rlmb"))

prior_list_normal <- lapply(ps_multi, function(ps) {
  list(mu = as.numeric(ps$mu), Sigma = ps$Sigma, dispersion = ps$dispersion)
})
out_normal <- do.call(multi_rNormal_reg, c(common, list(prior_list = prior_list_normal)))
stopifnot(inherits(out_normal, "mrglmb"), inherits(out_normal[[1L]], "rglmb"))

prior_list_ng <- lapply(ps_multi, function(ps) {
  list(
    mu = as.numeric(ps$mu), Sigma = ps$Sigma_0,
    shape = ps$shape, rate = ps$rate
  )
})
out_ng <- do.call(multi_rNormalGamma_reg, c(common, list(prior_list = prior_list_ng)))
stopifnot(inherits(out_ng, "mrglmb"))

prior_list_ing <- lapply(ps_multi, function(ps) {
  list(
    mu = as.numeric(ps$mu), Sigma = ps$Sigma,
    shape = ps$shape_ING, rate = ps$rate
  )
})
out_ing <- do.call(multi_rindepNormalGamma_reg, c(common, list(prior_list = prior_list_ing)))
stopifnot(inherits(out_ing, "mrglmb"))

for (out in list(out_rlmb, out_normal, out_ng, out_ing)) {
  for (j in seq_len(ncol(y))) {
    stopifnot(nrow(out[[j]]$coefficients) == n_draw, ncol(out[[j]]$coefficients) == 3L)
  }
}
cat("multi_* (iris): OK\n")

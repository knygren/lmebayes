# Smoke test: multi-response lmb on iris (glmbayes::lmb with cbind LHS)

if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install pkgload.")
pkgload::load_all(export_all = FALSE)

data("iris", package = "datasets")

set.seed(42)
n_draw <- 50L

form_multi <- cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species

ps_multi <- multi_prior_setup(form_multi, data = iris, family = gaussian())
stopifnot(inherits(ps_multi, "multi_PriorSetup"), length(ps_multi) == 4L)

pfamily_list <- lapply(ps_multi, function(ps) {
  dNormal_Gamma(
    mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate
  )
})

out <- lmb(
  form_multi,
  pfamily = pfamily_list,
  data = iris,
  n = n_draw,
  use_parallel = FALSE
)

stopifnot(inherits(out, "mlmb"), length(out) == 4L)
stopifnot(inherits(out[[1L]], "lmb"))
stopifnot(!inherits(out, "mrglmb"))
stopifnot(nrow(out[[1L]]$coefficients) == n_draw)

cm <- do.call(cbind, lapply(out, `[[`, "coef.means"))
stopifnot(nrow(cm) == 3L, ncol(cm) == 4L)

print(out)
s <- summary(out)
stopifnot(inherits(s, "summary.mlmb"), length(s) == 4L)
print(s)

cat("multi-response lmb (iris): OK\n")

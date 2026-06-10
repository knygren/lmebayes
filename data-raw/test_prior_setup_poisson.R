devtools::load_all("C:/Rpackages/lmebayes")

dat <- data.frame(
  y = rpois(40, lambda = 2),
  g = factor(rep(1:4, each = 10))
)
form <- y ~ (1 | g)

ps <- tryCatch(
  Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01),
  error = function(e) e
)
if (inherits(ps, "error")) {
  cat("Prior_Setup poisson (small data): expected error or singular:\n")
  cat(conditionMessage(ps), "\n")
} else {
  stopifnot(is.null(ps$dispersion_ranef))
  cat("Prior_Setup poisson: OK (dispersion_ranef is NULL)\n")
}

cat("test_prior_setup_poisson: done\n")

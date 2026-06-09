## Smoke test for multi_rNormal_reg
## Run from the package root: Rscript data-raw/test_multi_rNormal_reg.R

devtools::load_all(".", quiet = TRUE)

cat("---  multi_rNormal_reg  ---\n")
cat("formals:", paste(names(formals(multi_rNormal_reg)), collapse = ", "), "\n\n")

set.seed(42)
n_obs <- 50

## ---- input validation: list-x path ----

Y  <- cbind(y1 = rnorm(n_obs), y2 = rnorm(n_obs))
x1 <- cbind(1, rnorm(n_obs))               # p = 2
x2 <- cbind(1, rnorm(n_obs), rnorm(n_obs)) # p = 3
pl1 <- list(mu = c(0, 0),    Sigma = diag(2), dispersion = 1)
pl2 <- list(mu = c(0, 0, 0), Sigma = diag(3), dispersion = 1)

tryCatch(
  multi_rNormal_reg(1, Y, list(x1), list(pl1, pl2)),
  error = function(e) cat("[PASS] length mismatch caught:", conditionMessage(e), "\n")
)

tryCatch(
  multi_rNormal_reg(1, Y, list(x1, x2), pl1),
  error = function(e) cat("[PASS] single-prior misuse caught:", conditionMessage(e), "\n")
)

tryCatch(
  multi_rNormal_reg(1, Y, list(x1, x2), list(pl1)),
  error = function(e) cat("[PASS] prior_list wrong length caught:", conditionMessage(e), "\n")
)

## ---- list-x sampling ----

cat("\nSampling with list-x (n=1, varying p) ...\n")
res <- multi_rNormal_reg(
  n = 1, y = Y, x = list(x1, x2), prior_list = list(pl1, pl2),
  progbar = FALSE
)
stopifnot(is.list(res), !inherits(res, "mrglmb"))
stopifnot(identical(names(res), c("y1", "y2")))
stopifnot(ncol(res[["y1"]]$coefficients) == 2L)
stopifnot(ncol(res[["y2"]]$coefficients) == 3L)
cat("[PASS] list-x: class =", class(res),
    "| y1 coef dim =", dim(res[["y1"]]$coefficients),
    "| y2 coef dim =", dim(res[["y2"]]$coefficients), "\n")

## ---- shared-x path returns mrglmb ----

cat("\nSampling with shared-x (n=5, same p) ...\n")
x_s  <- cbind(1, rnorm(n_obs))
pl_s <- list(mu = c(0, 0), Sigma = diag(2), dispersion = 1)
res2 <- multi_rNormal_reg(
  n = 5, y = Y, x = x_s, prior_list = list(pl_s, pl_s),
  progbar = FALSE
)
stopifnot(inherits(res2, "mrglmb"))
stopifnot(all(dim(res2$coefficients) == c(5L, 2L)))
cat("[PASS] shared-x: class =", class(res2),
    "| coef dim =", dim(res2$coefficients), "\n")

cat("\nAll checks passed.\n")

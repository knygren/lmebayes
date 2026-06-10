# Smoke test: glmerb (Poisson) on bayesrules::airbnb with level-2 covariates
# and cross-level RE moderation (walk_c:rating_c, transit_c:log_price_c).
#
#   Rscript data-raw/test_glmerb_airbnb.R
#   Rscript data-raw/test_glmerb_airbnb.R quick
#   Rscript data-raw/test_glmerb_airbnb.R small

args <- commandArgs(trailingOnly = TRUE)
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))
use_small <- any(tolower(args) %in% c("small", "--small"))

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

if (use_small) {
  data("airbnb_small", package = "bayesrules")
  dat <- airbnb_small
  message("Using bayesrules::airbnb_small (n = ", nrow(dat), ")")
} else {
  data("airbnb", package = "bayesrules")
  dat <- airbnb
  message("Using bayesrules::airbnb (n = ", nrow(dat), ")")
}

dat$rating_c    <- dat$rating - mean(dat$rating)
dat$log_price_c <- scale(log(dat$price + 1))[, 1]
dat$walk_c      <- dat$walk_score - mean(dat$walk_score)
dat$transit_c   <- dat$transit_score - mean(dat$transit_score)
dat <- dat[complete.cases(dat[, c(
  "reviews", "rating", "rating_c", "price",
  "walk_score", "transit_score", "walk_c", "transit_c", "neighborhood"
)]), ]

form <- reviews ~
  walk_c + transit_c +
  rating_c + log_price_c +
  walk_c:rating_c + transit_c:log_price_c +
  (1 + rating_c + log_price_c || neighborhood)

design <- model_setup(form, data = dat, family = poisson())
stopifnot(isTRUE(design$rank_ok))
stopifnot(length(design$re_coef_names) == 3L)
stopifnot(ncol(design$X_hyper[["(Intercept)"]]) >= 3L)
stopifnot(ncol(design$X_hyper[["rating_c"]]) >= 2L)
stopifnot(ncol(design$X_hyper[["log_price_c"]]) >= 2L)

ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)
stopifnot(inherits(ps, "lmebayes_prior_setup"))
stopifnot(is.null(ps$dispersion_ranef))
stopifnot(all(diag(ps$Sigma_ranef) > 0))

n_draw <- if (run_quick) 50L else 200L
message("Posterior draws: n = ", n_draw)

fit <- glmerb(
  form,
  data = dat,
  family = poisson(),
  measurement_prior_list = ps,
  n = n_draw,
  seed = 42L
)

stopifnot(inherits(fit, "glmerb"))
stopifnot(length(fit$coef.mode) == 3L)
re_names <- fit$model_setup$re_coef_names
stopifnot(nrow(fit$fixef_draws[[re_names[1L]]]) == n_draw)

print(summary(fit))

cat("\ntest_glmerb_airbnb: OK\n")

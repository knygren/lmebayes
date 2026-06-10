# Smoke test: glmerb (Poisson) on bayesrules::airbnb neighborhood RE model.
#
# reviews ~ walkability covariates (level-2) + (1 | neighborhood).
# Same data prep as inst/examples/Ex_glmerb.R.
#
#   Rscript data-raw/test_glmerb_airbnb.R
#   Rscript data-raw/test_glmerb_airbnb.R quick    # n = 50
#   Rscript data-raw/test_glmerb_airbnb.R small    # airbnb_small

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

dat <- dat[complete.cases(dat[, c(
  "reviews", "neighborhood", "walk_score", "transit_score", "bike_score"
)]), ]
dat$walk_c    <- dat$walk_score    - mean(dat$walk_score)
dat$transit_c <- dat$transit_score - mean(dat$transit_score)
dat$bike_c    <- dat$bike_score    - mean(dat$bike_score)

form <- reviews ~ walk_c + transit_c + bike_c + (1 | neighborhood)

design <- model_setup(form, data = dat, family = poisson())
stopifnot(isTRUE(design$rank_ok))

ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)
stopifnot(inherits(ps, "lmebayes_prior_setup"))
stopifnot(is.null(ps$dispersion_ranef))

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
stopifnot(inherits(fit$glmer, "glmerMod"))
stopifnot(is.null(fit$lmer))
stopifnot(identical(fit$family$family, "poisson"))
re_names <- fit$model_setup$re_coef_names
stopifnot(nrow(fit$fixef_draws[[re_names[1L]]]) == n_draw)
stopifnot(nrow(fit$coefficients) == n_draw * nlevels(design$groups))
stopifnot(length(fit$coef.mode) == length(ps$prior_list))

print(summary(fit))

cat("\ntest_glmerb_airbnb: OK\n")

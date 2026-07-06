# Manual validation: glmerb Block~2 fixef + RE (glmer_full vs chain mean).
# Run after major sampler / pilot / ordering changes.
#
# Sampler: glmerb() defaults (tv_tol, gap_tol, mode_gap_max). No set.seed().
# Main draws: n = 3000 (full); n = 1000 with trailing arg "small" (both sections).
# progbar = TRUE only.
#
# Routes:
#   rGLMM_reg_known_vcov        (Poisson airbnb; Ex_13/14-style)
#   rGLMM_reg_estimated_vcov    (same Poisson airbnb; ING Block~2)
#
#   Rscript tests/manual/test_glmerb_mer_re_validation.R
#   Rscript tests/manual/test_glmerb_mer_re_validation.R small

args <- commandArgs(trailingOnly = TRUE)
use_small <- any(tolower(args) %in% c("small", "--small"))

source("tests/manual/_load.R")
source("tests/manual/_glmerb_re_validate.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_MIN  <- 1000L
N_FULL <- 3000L

## --- 1. Poisson known vcov (airbnb) ------------------------------------------

if (use_small) {
  data("airbnb_small", package = "bayesrules")
  dat <- airbnb_small
  n_poisson <- N_MIN
  message("Poisson: bayesrules::airbnb_small; main draws n = ", n_poisson)
} else {
  data("airbnb", package = "bayesrules")
  dat <- airbnb
  n_poisson <- N_FULL
  message("Poisson: bayesrules::airbnb; main draws n = ", n_poisson)
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

ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)

cat("\n=== glmerb Poisson known vcov; n =", n_poisson, "===\n\n")
fit <- glmerb(
  form,
  data         = dat,
  family       = poisson(),
  pfamily_list = pfamily_list(ps),
  n            = n_poisson,
  progbar      = TRUE
)

stopifnot(inherits(fit, "glmerb"))
re_names <- fit$model_setup$re_coef_names
stopifnot(identical(nrow(fit$fixef[[re_names[1L]]]), n_poisson))

lmebayes:::.validate_lmerb_block2_fixef_lmer(fit, label = "Poisson airbnb")

.validate_glmerb_re(fit, label = "Poisson airbnb")

## --- 2. ING estimated vcov (same Poisson airbnb model) -----------------------

ps_ing <- Prior_Setup_lmebayes(
  form,
  data           = dat,
  family         = poisson(),
  pwt            = 0.01
)
pf_ing <- pfamily_list(ps_ing, ptypes = "dIndependent_Normal_Gamma")

cat("\n=== glmerb ING estimated vcov (Poisson airbnb); n =", n_poisson, "===\n\n")
fit2 <- glmerb(
  form,
  data         = dat,
  family       = poisson(),
  pfamily_list = pf_ing,
  n            = n_poisson,
  progbar      = TRUE
)

stopifnot(isTRUE(fit2$prior$any_non_normal))
stopifnot(fit2$pilot_chisq$n_pilot > 0L)
stopifnot(is.matrix(fit2$fixef.dispersion))
stopifnot(all(is.finite(fit2$fixef.dispersion)), all(fit2$fixef.dispersion > 0))
stopifnot(identical(
  nrow(fit2$fixef[[fit2$model_setup$re_coef_names[1L]]]),
  n_poisson
))
re_names2 <- fit2$model_setup$re_coef_names
stopifnot(identical(re_names2, re_names))

lmebayes:::.validate_lmerb_block2_fixef_lmer(fit2, label = "Poisson airbnb ING")

.validate_glmerb_re(
  fit2,
  label = "Poisson airbnb ING",
  cor_icm_min = 0.85,
  cor_full_min = 0.8
)

cat("\ntest_glmerb_mer_re_validation: OK\n")

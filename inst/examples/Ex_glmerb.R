## glmerb: neighborhood random effects on bayesrules::airbnb_small (Poisson)
##
## Review counts for Chicago Airbnb listings.  Listing rating has a fixed
## main effect and a neighborhood random slope; neighborhood walkability
## (walk_c, constant within neighborhood) is a level-2 predictor of the
## neighborhood intercepts.  Poisson has no observation-level dispersion,
## so dispersion_ranef stays NULL.
##
## The full three-component workflow on the larger airbnb data (cross-level
## moderation, MCMC-vs-ICM diagnostics) is preserved as a demo:
##   demo("Ex_13_glmerb_Airbnb", package = "lmebayes")

if (requireNamespace("bayesrules", quietly = TRUE)) {

  data(airbnb_small, package = "bayesrules")
  dat <- airbnb_small
  dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
  dat$walk_c   <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
  dat <- dat[complete.cases(dat[, c("reviews", "rating_c", "walk_c",
                                    "neighborhood")]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))

  form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

  ## Default hyperpriors calibrated from a reference glmer fit (weak prior).
  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(),
                             pwt = 0.01)

  fit <- glmerb(
    form,
    data         = dat,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    n            = 1000L,
    seed         = 42L
  )
  summary(fit)

  ## Level-2 posterior means alongside the classical glmer reference.
  fit$coef.means
  lme4::fixef(fit$glmer)
}

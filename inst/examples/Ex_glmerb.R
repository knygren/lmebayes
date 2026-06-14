## glmerb: neighborhood random effects on bayesrules::airbnb_small (Poisson)
##
## Fast man-page example: prior setup, glmer reference fit, and ICM posterior
## mode only (simulate = FALSE).  No Gibbs draws — suitable for R CMD check.
##
## The full sampling workflow (pilot stage, main draws, centering diagnostics)
## is preserved as a demo:
##   demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")

if (requireNamespace("bayesrules", quietly = TRUE)) {

  data(airbnb_small, package = "bayesrules")
  dat <- airbnb_small
  dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
  dat$walk_c   <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
  dat <- dat[complete.cases(dat[, c("reviews", "rating_c", "walk_c",
                                    "neighborhood")]), ]
  dat$neighborhood <- droplevels(factor(dat$neighborhood))

  form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)

  fit <- glmerb(
    form,
    data         = dat,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    simulate     = FALSE
  )

  print_coef_means(fit)
  print(fit)
  lme4::fixef(fit$glmer)
}

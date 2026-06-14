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
  ps <- Prior_Setup_lmebayes(form, data = dat, family = poisson(), pwt = 0.01)
  
  # Inflate tau^2 by e.g. 25x to stress-test
  ps$Sigma_ranef <- ps$Sigma_ranef * 25
  

  fit <- glmerb(
    form,
    data         = dat,
    family       = poisson(),
    pfamily_list = pfamily_list(ps),
    n            = 2000L,
    seed         = 42L
  )
  cat("m_convergence used:", fit$convergence$m_convergence, "\n")
  cat(sprintf(
    "Pilot vs mode (chi-squared): p = %.4g\n",
    fit$pilot_mode_test$p_value
  ))

  ## Overall centering diagnostics: posterior mean vs pilot mean and mode.
  re_names <- fit$model_setup$re_coef_names
  X <- do.call(cbind, lapply(re_names, function(k) fit$fixef_draws[[k]]))
  cn <- unlist(lapply(re_names, function(k) {
    paste0(k, "::", colnames(fit$fixef_draws[[k]]))
  }))
  colnames(X) <- cn
  beta_bar <- colMeans(X)
  theta_pilot <- unlist(lapply(re_names, function(k) fit$coef.pilot.mean[[k]]))
  theta_mode  <- unlist(lapply(re_names, function(k) fit$coef.mode[[k]]))
  names(theta_pilot) <- cn
  names(theta_mode)  <- cn

  S <- stats::cov(X)
  V <- S / nrow(X)
  V_inv <- solve(V)
  d_pilot <- beta_bar - theta_pilot
  d_mode  <- beta_bar - theta_mode
  Q_pilot <- as.numeric(t(d_pilot) %*% V_inv %*% d_pilot)
  Q_mode  <- as.numeric(t(d_mode)  %*% V_inv %*% d_mode)
  p_pilot <- stats::pchisq(Q_pilot, df = ncol(X), lower.tail = FALSE)
  p_mode  <- stats::pchisq(Q_mode,  df = ncol(X), lower.tail = FALSE)
  cat(sprintf(
    "Overall centering (chi-squared): p(mean=pilot)=%.4g, p(mean=mode)=%.4g\n",
    p_pilot, p_mode
  ))

  ## Level-2 posterior means alongside the classical glmer reference.
  print_coef_means(fit)
  lme4::fixef(fit$glmer)

  ## Full print and summary.
  print(fit)
  summary(fit)
}

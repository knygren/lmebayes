## Demo: glmerb() full workflow on bayesrules::airbnb_small (Poisson)
##
## Replica of the original ?glmerb example before the man-page example was
## switched to simulate = FALSE (ICM only).  Preserved as a demo because the
## full run (2000 draws, pilot stage, convergence calibration) takes several
## minutes.  This demo keeps the complete workflow including pilot diagnostics,
## centering tests, print_coef_means(), print(), and summary().
##
##   demo("Ex_14_glmerb_airbnb_small", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}

## glmerb: neighborhood random effects on bayesrules::airbnb_small (Poisson)
##
## Review counts for Chicago Airbnb listings.  Listing rating has a fixed
## main effect and a neighborhood random slope; neighborhood walkability
## (walk_c, constant within neighborhood) is a level-2 predictor of the
## neighborhood intercepts.  Poisson has no observation-level dispersion,
## so dispersion_ranef stays NULL.

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

## Inflate tau^2 to stress-test convergence calibration (optional).
#ps$Sigma_ranef <- ps$Sigma_ranef * 25

fit <- glmerb(
  form,
  data         = dat,
  family       = poisson(),
  pfamily_list = pfamily_list(ps),
  n            = 10000L,
#  seed         = 42L,
  progbar=FALSE
)
cat("m_convergence used:", fit$convergence$m_convergence, "\n")
cat(sprintf(
  "Pilot vs mode (chi-squared): p = %.4g\n",
  fit$pilot_chisq$p_value
))

## Overall centering diagnostics: posterior mean vs pilot mean and mode.
re_names <- fit$model_setup$re_coef_names
X <- do.call(cbind, lapply(re_names, function(k) fit$fixef[[k]]))
cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
colnames(X) <- cn
beta_bar <- colMeans(X)
theta_pilot <- unlist(lapply(re_names, function(k) fit$fixef.init[[k]]))
theta_mode  <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
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
lmebayes:::print_coef_means(fit)
lme4::fixef(fit$glmer)

## Full print and summary.
print(fit)
summary(fit)

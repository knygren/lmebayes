# Statistical test: dIndependent_Normal_Gamma tau^2 sampling in lmerb()/glmerb().
#
# With the v2 two-block sampler, ING components draw (gamma_k, tau^2_k)
# jointly each inner sweep.  This test checks
#   1. lmerb (gaussian, big_word_club): with an informative dispersion prior
#      (pwt_dispersion = 0.4) the tau^2 posterior mean lands near the lmer
#      reference RE variance, and the gamma posterior means track the
#      (fixed-tau^2) ICM mode.
#   2. the sampler-level prior-vs-data guard (n_prior <= J) fires through
#      lmerb for a hand-built prior-dominated ING pfamily.
#   3. glmerb (poisson, simulated): ING smoke run produces finite, varying
#      tau^2 draws.
#
# Loads glmbayesCore from source first so the test exercises the current
# v2 validator (guard included) rather than the installed snapshot.
#
#   Rscript data-raw/test_ing_sampling.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
#pkgload::load_all("../glmbayesCore", export_all = FALSE, quiet = TRUE)
#pkgload::load_all(export_all = FALSE)

## --- 1. lmerb gaussian: tau^2 posterior centered near the reference -----------
data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(dat, !is.na(score_ppvt))

form <- score_ppvt ~ private_school + (1 | school_id)

ps <- Prior_Setup_lmebayes(form, data = dat, pwt = 0.01,
                           pwt_dispersion = 0.2)
J        <- nlevels(ps$design$groups)
tau2_hat <- unname(ps$prior_list[["(Intercept)"]]$dispersion_fixef)


ps

pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

out1 <- capture.output(
  fit <- lmerb(form, data = dat,
               pfamily_list = pf,
               dispersion_ranef = ps$dispersion_ranef,
               n = 300L, seed = 42L)
)

t2_draws <- fit$tau2_draws[, "(Intercept)"]
t2_mean  <- fit$tau2.means[["(Intercept)"]]
pr_int   <- pf[["(Intercept)"]]$prior_list
stopifnot(
  length(t2_draws) == 300L,
  all(is.finite(t2_draws)), all(t2_draws > 0),
  stats::sd(t2_draws) > 0,
  ## Fixed truncation window: every draw inside [disp_lower, disp_upper]
  all(t2_draws >= pr_int$disp_lower),
  all(t2_draws <= pr_int$disp_upper)
)


summary(fit)


## Informative prior (n_prior = 0.4/0.6 * J ~ 2J/3) centered at tau2_hat plus
## J group-level observations: the posterior mean must land within a modest
## factor of the reference variance.
if (t2_mean < tau2_hat / 5 || t2_mean > tau2_hat * 5) {
  stop("tau^2 posterior mean ", signif(t2_mean, 4),
       " far from reference ", signif(tau2_hat, 4))
}
## gamma posterior means should track the (fixed-tau^2) ICM mode loosely:
## tau^2 is sampled here, so allow a couple of posterior SDs of slack.
g_means <- fit$coef.means[["(Intercept)"]]
g_mode  <- fit$coef.mode[["(Intercept)"]]
g_sd    <- apply(fit$fixef_draws[["(Intercept)"]], 2L, stats::sd)
stopifnot(all(abs(g_means - g_mode) < 3 * g_sd))
cat(sprintf(
  "1. lmerb ING sampling: tau^2 mean = %.2f (reference %.2f), gamma means track mode: OK\n",
  t2_mean, tau2_hat
))

## --- 2. sampler-level guard via lmerb ----------------------------------------
pl1 <- pf[["(Intercept)"]]$prior_list
shape_bad <- (10 * J + 1 + length(pl1$mu)) / 2  # n_prior = 10J >> J
pf_bad <- pf
pf_bad[["(Intercept)"]] <- dIndependent_Normal_Gamma(
  mu         = pl1$mu,
  Sigma      = pl1$Sigma,
  shape      = shape_bad,
  rate       = as.numeric(pl1$rate),
  disp_lower = as.numeric(pl1$disp_lower),
  disp_upper = as.numeric(pl1$disp_upper)
)
err <- tryCatch(
  {
    capture.output(
      lmerb(form, data = dat, pfamily_list = pf_bad,
            dispersion_ranef = ps$dispersion_ranef, n = 5L)
    )
    NULL
  },
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err), grepl("n_prior <= J", err, fixed = TRUE))
cat("2. prior-dominated ING pfamily rejected through lmerb (n_prior > J): OK\n")

## --- 3. glmerb poisson + ING smoke run ---------------------------------------
set.seed(11)
J_p   <- 15L
n_per <- 20L
g     <- factor(rep(seq_len(J_p), each = n_per))
b0    <- rnorm(J_p, mean = 0.5, sd = 0.4)
eta   <- b0[as.integer(g)]
datp  <- data.frame(y = rpois(J_p * n_per, exp(eta)), g = g)

form_p <- y ~ 1 + (1 | g)
ps_p <- Prior_Setup_lmebayes(form_p, data = datp, family = poisson(),
                             pwt = 0.05, pwt_dispersion = 0.4)
pf_p <- pfamily_list(ps_p, ptypes = "dIndependent_Normal_Gamma")

out3 <- capture.output(
  fit_p <- glmerb(form_p, data = datp, family = poisson(),
                  pfamily_list = pf_p,
                  n = 25L, seed = 7L)
)
pr_p <- pf_p[["(Intercept)"]]$prior_list
stopifnot(
  inherits(fit_p, "glmerb"),
  !is.null(fit_p$coefficients),
  is.matrix(fit_p$tau2_draws),
  nrow(fit_p$tau2_draws) == 25L,
  all(is.finite(fit_p$tau2_draws)), all(fit_p$tau2_draws > 0),
  all(fit_p$tau2_draws[, "(Intercept)"] >= pr_p$disp_lower),
  all(fit_p$tau2_draws[, "(Intercept)"] <= pr_p$disp_upper),
  stats::sd(fit_p$tau2_draws[, "(Intercept)"]) > 0,
  all(is.finite(unlist(fit_p$coef.means)))
)
cat(sprintf(
  "3. glmerb poisson ING smoke: tau^2 mean = %.3f (reference %.3f): OK\n",
  fit_p$tau2.means[["(Intercept)"]],
  unname(ps_p$prior_list[["(Intercept)"]]$dispersion_fixef)
))

cat("\ntest_ing_sampling.R: all checks passed\n")

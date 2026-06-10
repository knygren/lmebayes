devtools::load_all("C:/Rpackages/lmebayes")

dat <- data.frame(
  y = rnorm(40),
  g = factor(rep(1:4, each = 10))
)
form <- y ~ (1 | g)

design <- model_setup(form, data = dat, family = gaussian())
stopifnot(inherits(design$lmer_fit, "lmerMod"))

dat_pois <- dat
dat_pois$y <- rpois(nrow(dat), lambda = 2)
design2 <- model_setup(form, data = dat_pois, family = poisson())
stopifnot(inherits(design2$glmer_fit, "glmerMod"))

ps <- structure(
  list(
    formula = form,
    prior_list = list(
      `(Intercept)` = list(
        mu_fixef = c(`(Intercept)` = 0),
        Sigma_fixef = matrix(1),
        dispersion_fixef = 1
      )
    ),
    Sigma_ranef = diag(1),
    dispersion_ranef = 1
  ),
  class = "lmebayes_prior_setup"
)

fit <- glmerb(form, data = dat, family = gaussian(),
              measurement_prior_list = ps, n = 2, simulate = FALSE)
stopifnot(inherits(fit$glmer, c("glmerMod", "lmerMod")))
stopifnot(is.null(fit$lmer))
stopifnot(!is.null(ps$dispersion_ranef))

bad <- tryCatch(
  model_setup(y ~ x + (1 + x | g), data = dat, family = gaussian()),
  error = function(e) e
)
stopifnot(inherits(bad, "error"))

bad2 <- tryCatch(
  model_setup(y ~ (1 || g), data = dat, family = gaussian()),
  error = function(e) e
)
stopifnot(inherits(bad2, "error"))

cat("test_glmerb_model_setup: OK\n")

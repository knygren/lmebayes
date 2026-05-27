#!/usr/bin/env Rscript

## Gamma–Gamma conjugate: Gamma prior on the *rate* β of a Gamma likelihood
## (intercept-only, identity link on β).
##
## **Model**
##   Y_i | β ~ Gamma(shape = k, rate = β),  i = 1,...,n  (k known)
##   β ~ Gamma(α₀, β₀)                                   (prior)
##   β | y ~ Gamma(α₀ + n·k, β₀ + Σyᵢ)                  (posterior)
##
## With k = 1 the Gamma response collapses to the Exponential distribution,
## giving the simplest possible closed-form update:
##   β | y ~ Gamma(α₀ + n, β₀ + Σyᵢ).
##
## **Interpretation of the "identity" link here**
##   dGamma(Inv_Dispersion=FALSE) with family = Gamma(link = "identity") treats the
##   intercept coefficient as the *rate* β directly (not the mean μ = k/β).
##   This matches the conjugate table: "Gamma with known shape | rate →
##   Gamma conjugate prior on rate".
##
## **Examples**
##   A. Simulated exponential data (k = 1): exact analytic check.
##   B. cherry_blossom bloom-day data from bayesrules (k estimated from data,
##      used as a fixed shape to illustrate the workflow with real data).
##
## Suggested packages: bayesrules, ggplot2.
## Run after installing the package, e.g.:
##   Rscript inst/scripts/bayesrules_gamma_gamma_glmb.R

library(glmbayes)

if (requireNamespace("bayesrules", quietly = TRUE)) {
  library(bayesrules)
  has_br <- TRUE
} else {
  message("bayesrules not installed; Example B will be skipped.")
  has_br <- FALSE
}
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  has_gg <- TRUE
} else {
  has_gg <- FALSE
}


## =============================================================================
## Example A — simulated exponential data  (k = 1, closed-form is exact)
## =============================================================================

cat("\n=== Example A: Exponential data (Gamma shape k = 1) ===\n\n")

## True rate β = 0.5  →  mean waiting time = 1/β = 2
set.seed(42)
n_A   <- 30L
beta_true <- 0.5
y_A   <- rexp(n_A, rate = beta_true)

## Prior: β ~ Gamma(α₀ = 2, β₀ = 4)  →  prior mean rate = 2/4 = 0.5
alpha0_A <- 2
beta0_A  <- 4

## Analytic conjugate posterior (k = 1):
##   shape_post = α₀ + n·k = 2 + 30 = 32
##   rate_post  = β₀ + Σyᵢ = 4 + sum(y_A)
shape_post_A <- alpha0_A + n_A * 1
rate_post_A  <- beta0_A  + sum(y_A)

cat(sprintf("Observed data:  n = %d,  mean(y) = %.4f,  sum(y) = %.4f\n",
            n_A, mean(y_A), sum(y_A)))
cat(sprintf("Prior:          Gamma(shape = %g, rate = %g)  →  prior mean rate = %.4f\n",
            alpha0_A, beta0_A, alpha0_A / beta0_A))
cat(sprintf("Posterior:      Gamma(shape = %.1f, rate = %.4f)\n",
            shape_post_A, rate_post_A))
cat(sprintf("  Posterior mean rate  = %.6f  (true β = %.2f)\n",
            shape_post_A / rate_post_A, beta_true))
cat(sprintf("  Posterior mean μ(=1/β) = %.6f  (true mean = %.2f)\n",
            rate_post_A / shape_post_A, 1 / beta_true))
cat(sprintf("  90%% credible interval for β: [%.4f, %.4f]\n\n",
            qgamma(0.05, shape_post_A, rate_post_A),
            qgamma(0.95, shape_post_A, rate_post_A)))

## Set up the pfamily — lik_shape = 1 (exponential / Gamma with known shape 1)
beta_init_A <- matrix(alpha0_A / beta0_A, nrow = 1L, ncol = 1L)
colnames(beta_init_A) <- "(Intercept)"

pf_A <- dGamma(
  shape          = alpha0_A,
  rate           = beta0_A,
  beta           = beta_init_A,
  Inv_Dispersion = FALSE,
  lik_shape      = 1          ## k = 1  (exponential)
)

## Fit with glmb() — family = Gamma(link = "identity"), coefficient = rate β
data_A <- data.frame(y = y_A)
set.seed(2026)
fit_A <- glmb(
  n       = 20000,
  y ~ 1,
  data    = data_A,
  family  = Gamma(link = "identity"),
  pfamily = pf_A
)

cat("glmb() summary (coefficient = Gamma rate β):\n")
print(summary(fit_A))

## Verify: glmb draws should match the analytic posterior
smp_A <- fit_A$coefficients[, 1L]
cat(sprintf("glmb draw mean = %.6f  |  analytic posterior mean = %.6f\n",
            mean(smp_A), shape_post_A / rate_post_A))
cat(sprintf("glmb draw SD   = %.6f  |  analytic posterior SD   = %.6f\n",
            sd(smp_A), sqrt(shape_post_A) / rate_post_A))

## Overlay draws on analytic density
if (has_gg) {
  grid_A <- seq(
    qgamma(0.001, shape_post_A, rate_post_A),
    qgamma(0.999, shape_post_A, rate_post_A),
    length.out = 400
  )
  p_A <- ggplot(data.frame(beta = smp_A), aes(beta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 55,
                   fill = "steelblue", alpha = 0.4, color = NA) +
    geom_line(
      data = data.frame(beta = grid_A,
                        d    = dgamma(grid_A, shape_post_A, rate_post_A)),
      aes(beta, d), linewidth = 1, colour = "black", linetype = "dashed",
      inherit.aes = FALSE
    ) +
    labs(
      title    = "Gamma–Gamma conjugate: glmb draws vs analytic posterior",
      subtitle = sprintf(
        "Exponential data (k=1, true β=%.2f); prior Gamma(%g,%g); n=%d",
        beta_true, alpha0_A, beta0_A, n_A
      ),
      x = expression(beta ~ "(rate)"),
      y = "density"
    ) +
    theme_bw()
  print(p_A)
}


## =============================================================================
## Example B — cherry_blossom bloom days from bayesrules  (k estimated)
## =============================================================================

if (has_br) {
  cat("\n=== Example B: cherry_blossom bloom days (bayesrules) ===\n\n")

  ## cherry_blossom has `day` = Julian day of first bloom (~85–125)
  data(cherry_blossom, package = "bayesrules")
  y_B <- cherry_blossom$day
  y_B <- y_B[!is.na(y_B) & y_B > 0]   ## keep positive non-missing values
  n_B <- length(y_B)

  cat(sprintf("cherry_blossom: n = %d,  mean = %.2f,  SD = %.2f\n",
              n_B, mean(y_B), sd(y_B)))

  ## Estimate Gamma shape from method-of-moments: k ≈ mean²/var
  k_est_B <- mean(y_B)^2 / var(y_B)
  cat(sprintf("Method-of-moments shape estimate: k = %.3f\n", k_est_B))
  cat("  (We fix k at this value for the conjugate update.)\n\n")

  ## Prior: β ~ Gamma(α₀, β₀) with prior mean rate = k/mean(y)
  ## Use a diffuse prior: α₀ = 2, β₀ chosen so prior mean rate = k_est / mean(y)
  alpha0_B <- 2
  beta0_B  <- alpha0_B / (k_est_B / mean(y_B))   ## prior mean rate = k/mean(y)

  ## Analytic posterior
  shape_post_B <- alpha0_B + n_B * k_est_B
  rate_post_B  <- beta0_B  + sum(y_B)

  cat(sprintf("Prior:     Gamma(shape = %.3f, rate = %.6f)  →  prior mean rate = %.6f\n",
              alpha0_B, beta0_B, alpha0_B / beta0_B))
  cat(sprintf("Posterior: Gamma(shape = %.3f, rate = %.4f)\n",
              shape_post_B, rate_post_B))
  cat(sprintf("  Posterior mean rate  = %.6f\n", shape_post_B / rate_post_B))
  cat(sprintf("  Posterior mean bloom day (k/β) = %.4f\n",
              k_est_B / (shape_post_B / rate_post_B)))
  cat(sprintf("  90%% credible interval for rate β: [%.6f, %.6f]\n\n",
              qgamma(0.05, shape_post_B, rate_post_B),
              qgamma(0.95, shape_post_B, rate_post_B)))

  beta_init_B <- matrix(alpha0_B / beta0_B, nrow = 1L, ncol = 1L)
  colnames(beta_init_B) <- "(Intercept)"

  pf_B <- dGamma(
    shape          = alpha0_B,
    rate           = beta0_B,
    beta           = beta_init_B,
    Inv_Dispersion = FALSE,
    lik_shape      = k_est_B    ## fixed Gamma shape from method-of-moments
  )

  data_B <- data.frame(y = y_B)
  set.seed(2026)
  fit_B <- glmb(
    n       = 20000,
    y ~ 1,
    data    = data_B,
    family  = Gamma(link = "identity"),
    pfamily = pf_B
  )

  cat("glmb() summary (cherry_blossom, coefficient = Gamma rate β):\n")
  print(summary(fit_B))

  smp_B <- fit_B$coefficients[, 1L]
  cat(sprintf("glmb draw mean = %.6f  |  analytic posterior mean = %.6f\n",
              mean(smp_B), shape_post_B / rate_post_B))

  if (has_gg) {
    grid_B <- seq(
      qgamma(0.001, shape_post_B, rate_post_B),
      qgamma(0.999, shape_post_B, rate_post_B),
      length.out = 400
    )
    p_B <- ggplot(data.frame(beta = smp_B), aes(beta)) +
      geom_histogram(aes(y = after_stat(density)), bins = 55,
                     fill = "darkorange", alpha = 0.4, color = NA) +
      geom_line(
        data = data.frame(beta = grid_B,
                          d    = dgamma(grid_B, shape_post_B, rate_post_B)),
        aes(beta, d), linewidth = 1, colour = "black", linetype = "dashed",
        inherit.aes = FALSE
      ) +
      labs(
        title    = "Gamma–Gamma conjugate: glmb draws vs analytic posterior",
        subtitle = sprintf(
          "cherry_blossom bloom days; fixed k = %.2f; n = %d",
          k_est_B, n_B
        ),
        x = expression(beta ~ "(rate)"),
        y = "density"
      ) +
      theme_bw()
    print(p_B)
  }
}

message("\nFinished bayesrules + glmbayes Gamma–Gamma conjugate demos.")

# Precompute Eight Schools two-block Gibbs outputs for vignettes/Chapter-13.Rmd:
#   (1) Population block: dNormal_Gamma (Gibbs-suitable)
#   (2) Population block: dIndependent_Normal_Gamma (slower)
# Same seeds, burn-in, and stored iterations as the vignette.
#
# Run from package root:
#   Rscript data-raw/make_Chapter13_Eight_Schools_gibbs_output.R
# Optional path:
#   Rscript data-raw/make_Chapter13_Eight_Schools_gibbs_output.R "C:/path/to/glmbayes"
#
# Writes: inst/extdata/Chapter13_Eight_Schools_two_gibbs.rds

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload, e.g. install.packages('pkgload')")
}

pkgload::load_all(export_all = FALSE)

## --- schools-data (same as vignette Chapter-13.Rmd) ---
school   <- c("A", "B", "C", "D", "E", "F", "G", "H")
estimate <- c(28.39, 7.94, -2.75, 6.82, -0.64, 0.63, 18.01, 12.16)
sd_obs   <- c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6)
J        <- length(school)
sigma_y_sq <- sd_obs^2

mu_mu    <- mean(estimate)
sigma_mu <- var(estimate)
n_prior  <- 0.5
disp_ML  <- var(estimate)
shape    <- n_prior / 2
rate     <- disp_ML * shape

x_one <- as.matrix(rep(1, J), nrow = J, ncol = 1)

## ========== (1) dNormal_Gamma ==========
message("Eight Schools: dNormal_Gamma Gibbs ...")
set.seed(101)
theta_ng <- estimate
n_burn_ng <- 1000L
n_sim_ng  <- 1000L
theta_out_ng <- matrix(0, nrow = n_sim_ng, ncol = J)
mu_out_ng <- numeric(n_sim_ng)
sigma_theta_out_ng <- numeric(n_sim_ng)

for (k in seq_len(n_burn_ng)) {
  out_pop <- rlmb(
    1, y = theta_ng, x = x_one,
    pfamily = dNormal_Gamma(mu_mu, sigma_mu / disp_ML, shape = shape, rate = rate)
  )
  mu_theta       <- out_pop$coefficients[1, 1]
  sigma_theta_sq <- out_pop$dispersion
  for (j in seq_len(J)) {
    theta_ng[j] <- rlmb(
      1, y = estimate[j], x = as.matrix(1),
      pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
    )$coefficients[1, 1]
  }
}

for (k in seq_len(n_sim_ng)) {
  out_pop <- rlmb(
    1, y = theta_ng, x = x_one,
    pfamily = dNormal_Gamma(mu_mu, sigma_mu / disp_ML, shape = shape, rate = rate)
  )
  mu_theta       <- out_pop$coefficients[1, 1]
  sigma_theta_sq <- out_pop$dispersion
  for (j in seq_len(J)) {
    theta_ng[j] <- rlmb(
      1, y = estimate[j], x = as.matrix(1),
      pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
    )$coefficients[1, 1]
  }
  theta_out_ng[k, ]     <- theta_ng
  mu_out_ng[k]          <- mu_theta
  sigma_theta_out_ng[k] <- sqrt(sigma_theta_sq)
}
colnames(theta_out_ng) <- school

## ========== (2) dIndependent_Normal_Gamma ==========
message("Eight Schools: dIndependent_Normal_Gamma Gibbs ...")
theta <- estimate
n_burn_schools <- 1000L
n_sim_schools  <- 1000L
theta_out <- matrix(0, nrow = n_sim_schools, ncol = J)
mu_out    <- numeric(n_sim_schools)
sigma_theta_out <- numeric(n_sim_schools)
iters_out1 <- numeric(n_burn_schools)
iters_out2 <- numeric(n_sim_schools)

set.seed(102)
for (k in seq_len(n_burn_schools)) {
  out_pop <- rlmb(
    1, y = theta, x = x_one,
    pfamily = dIndependent_Normal_Gamma(mu_mu, sigma_mu, shape = shape, rate = rate)
  )
  mu_theta       <- out_pop$coefficients[1, 1]
  sigma_theta_sq <- out_pop$dispersion
  for (j in seq_len(J)) {
    theta[j] <- rlmb(
      1, y = estimate[j], x = as.matrix(1),
      pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
    )$coefficients[1, 1]
  }
  iters_out1[k] <- out_pop$iters
}

for (k in seq_len(n_sim_schools)) {
  out_pop <- rlmb(
    1, y = theta, x = x_one,
    pfamily = dIndependent_Normal_Gamma(mu_mu, sigma_mu, shape = shape, rate = rate)
  )
  mu_theta       <- out_pop$coefficients[1, 1]
  sigma_theta_sq <- out_pop$dispersion
  for (j in seq_len(J)) {
    theta[j] <- rlmb(
      1, y = estimate[j], x = as.matrix(1),
      pfamily = dNormal(mu_theta, Sigma = sigma_theta_sq, dispersion = sigma_y_sq[j])
    )$coefficients[1, 1]
  }
  theta_out[k, ] <- theta
  mu_out[k]      <- mu_theta
  sigma_theta_out[k] <- sqrt(sigma_theta_sq)
  iters_out2[k]  <- out_pop$iters
}
colnames(theta_out) <- school

cmp_theta <- data.frame(
  school = school,
  raw = estimate,
  NG_mean = colMeans(theta_out_ng),
  NG_SD   = apply(theta_out_ng, 2, sd),
  Indep_mean = colMeans(theta_out),
  Indep_SD   = apply(theta_out, 2, sd)
)
cmp_hyper <- data.frame(
  parameter = c("mu", "sigma_theta", "mean_pop_iters_main"),
  NG_mean = c(mean(mu_out_ng), mean(sigma_theta_out_ng), NA_real_),
  Indep_mean = c(mean(mu_out), mean(sigma_theta_out), mean(iters_out2))
)

out <- list(
  normal_gamma = list(
    seed = 101L,
    n_burn = n_burn_ng,
    n_sim = n_sim_ng,
    theta_out = theta_out_ng,
    mu_out = mu_out_ng,
    sigma_theta_out = sigma_theta_out_ng
  ),
  indep_norm_gamma = list(
    seed = 102L,
    n_burn = n_burn_schools,
    n_sim = n_sim_schools,
    theta_out = theta_out,
    mu_out = mu_out,
    sigma_theta_out = sigma_theta_out,
    iters_out1 = iters_out1,
    iters_out2 = iters_out2,
    mean_iters_burn = mean(iters_out1),
    mean_iters_main = mean(iters_out2)
  ),
  comparison = list(by_school = cmp_theta, hyper = cmp_hyper),
  school = school,
  estimate = estimate,
  package_version = as.character(utils::packageVersion("glmbayes")),
  generated_at = Sys.time()
)

dir.create(file.path(root, "inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "Chapter13_Eight_Schools_two_gibbs.rds")
saveRDS(out, dest, compress = "xz")
message("Wrote ", normalizePath(dest, winslash = "/"))

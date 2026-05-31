# =============================================================================
# Self-contained neighborhood block Gibbs — Block 1: rindepNormalGamma_reg
# =============================================================================
#
# Same hierarchical model as benchmark_airbnb_neighborhood_rNormalGamma_reg_block.R
# but Block 1 uses the *independent* Normal–Gamma prior via
# rindepNormalGamma_reg_multi (matrix y, prior_list = list of per-column prior_lists).
#
#   Rscript data-raw/benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block.R
#   Rscript data-raw/benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block.R quick
#
# Companion (conjugate NG hyper): benchmark_airbnb_neighborhood_rNormalGamma_reg_block.R

# --- CLI ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_legacy <- any(tolower(args) %in% c("legacy", "--legacy", "-l"))
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))
use_opencl_fixed <- FALSE
use_opencl_random <- FALSE
if (any(tolower(args) %in% c("opencl_fixed", "use_opencl_fixed", "ocl_fixed"))) {
  use_opencl_fixed <- TRUE
}
if (any(tolower(args) %in% c("opencl_random", "use_opencl_random", "ocl_random"))) {
  use_opencl_random <- TRUE
}
path_args <- args[!tolower(args) %in% c(
  "legacy", "--legacy", "-l",
  "quick", "--quick", "-q",
  "opencl_fixed", "use_opencl_fixed", "ocl_fixed",
  "opencl_random", "use_opencl_random", "ocl_random"
)]

root <- if (length(path_args) >= 1L) {
  normalizePath(path_args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.")
}
pkgload::load_all(export_all = FALSE)

OUT_RDS <- "Airbnb_neighborhood_rindepNormalGamma_reg_benchmark.rds"

# --- Small helpers (this script only) ----------------------------------------
fmt_hms <- function(secs) {
  secs <- max(0, as.numeric(secs))
  h <- floor(secs / 3600)
  rem <- secs - h * 3600
  m <- floor(rem / 60)
  s <- rem - m * 60
  sprintf("%d h %d min %.2f s", h, m, s)
}

gibbs_report_interval <- function(n) {
  max(1L, min(50L, as.integer(n %/% 10L)))
}

gibbs_message_first_iter_estimates <- function(tag, sec_per_iter, n_burn, n_sim) {
  message(tag, "After iteration 1 — estimated burn-in: ",
          fmt_hms(sec_per_iter * n_burn),
          "; estimated total: ",
          fmt_hms(sec_per_iter * (n_burn + n_sim)),
          " (", signif(sec_per_iter, 4), " s/iter)")
}

.recycle_hyper_gamma <- function(x, l1, name) {
  x <- as.numeric(x)
  if (length(x) == 1L) return(rep(x, l1))
  if (length(x) == l1) return(x)
  stop("'", name, "' must have length 1 or ", l1, ".", call. = FALSE)
}

build_prior_list_vector <- function(ps, coef_names, max_disp_perc = 0.99) {
  l1 <- length(ps$mu)
  sig_diag <- diag(ps$Sigma)
  mu_num <- as.numeric(ps$mu)
  shape_v <- .recycle_hyper_gamma(ps$shape_ING, l1, "shape_ING")
  rate_v <- .recycle_hyper_gamma(ps$rate, l1, "rate")
  prior_lists <- lapply(seq_len(l1), function(j) {
    pl <- list(
      mu = c(mu_num[j]),
      Sigma = matrix(sig_diag[j], 1L, 1L),
      shape = shape_v[j],
      rate = rate_v[j],
      max_disp_perc = max_disp_perc
    )
    attr(pl, "Prior Type") <- "dIndependent_Normal_Gamma"
    pl
  })
  prior_lists
}

print_coda_hyper_chain <- function(hyper_out, col_names, label) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    message("\n--- CODA (", label, "): skipped — install.packages('coda') ---")
    return(NULL)
  }
  if (nrow(hyper_out) < 2L) {
    message("\n--- CODA (", label, "): skipped — need n_sim >= 2 ---")
    return(NULL)
  }
  mcmc_obj <- coda::mcmc(hyper_out)
  colnames(mcmc_obj) <- col_names
  message("\n--- CODA summary: ", label, " ---")
  print(summary(mcmc_obj))
  ess <- coda::effectiveSize(mcmc_obj)
  message("Effective sample size:")
  print(stats::setNames(ess, col_names))
  list(
    summary = summary(mcmc_obj),
    effective_size = as.list(ess)
  )
}

## Block 1 wrapper: population hyper (mu, Sigma) for Block 2 from ING vector draw.
block1_hyper_draw <- function(beta_mat,
                              prior_list_vector,
                              beta_names,
                              use_opencl_fixed = FALSE) {
  colnames(beta_mat) <- beta_names
  out <- rindepNormalGamma_reg_multi(
    n = 1L,
    y = beta_mat,
    x = matrix(1, k, 1, dimnames = list(NULL, "(Intercept)")),
    prior_list = prior_list_vector,
    family = gaussian(),
    use_parallel = FALSE,
    use_opencl = use_opencl_fixed,
    verbose = FALSE,
    progbar = FALSE
  )
  hyper_from_block1_draw(out, beta_names)
}

## Block 2 hyper: mu = intercept draw; Sigma = diag(ING dispersion draws) per column.
## out is an "mrglmb" object: named list of rglmb fits, one per column of y.
hyper_from_block1_draw <- function(out, beta_names) {
  mu   <- vapply(out, function(fit) fit$coefficients[1L, 1L], numeric(1))
  disp <- vapply(out, function(fit) fit$dispersion[1L],       numeric(1))
  Sigma <- diag(disp)
  dimnames(Sigma) <- list(beta_names, beta_names)
  list(
    mu = stats::setNames(mu, beta_names),
    Sigma = Sigma,
    dispersion = stats::setNames(disp, beta_names)
  )
}

# --- Data --------------------------------------------------------------------
n_burn <- if (run_quick) 5L else 200L
n_sim  <- if (run_quick) 10L else 1000L
n_gibbs <- n_burn + n_sim

if (run_quick) {
  data("airbnb_small", package = "bayesrules")
  airbnb_dat <- airbnb_small
  message("Quick mode: bayesrules::airbnb_small")
} else {
  data("airbnb", package = "bayesrules")
  airbnb_dat <- airbnb
}

airbnb_dat$rating_c <- airbnb_dat$rating - mean(airbnb_dat$rating)
airbnb_dat$room_type <- factor(airbnb_dat$room_type)
airbnb_dat <- airbnb_dat[complete.cases(airbnb_dat[, c(
  "reviews", "rating_c", "room_type", "neighborhood",
  "walk_score", "transit_score", "bike_score"
)]), ]

form_x <- reviews ~ rating_c + room_type
X_full <- model.matrix(form_x, data = airbnb_dat)
y_full <- airbnb_dat$reviews
block_full <- factor(airbnb_dat$neighborhood)
l2 <- length(y_full)
l1 <- ncol(X_full)
k <- nlevels(block_full)
beta_names <- colnames(X_full)
neighborhood_names <- levels(block_full)
block_info <- glmbayes:::normalize_block(block_full, l2)

## Neighborhood-level design (not used in Block 1 yet): intercept + centered scores
nbhd_unique <- airbnb_dat[!duplicated(airbnb_dat$neighborhood), c(
  "neighborhood", "walk_score", "transit_score", "bike_score"
)]
nbhd_unique <- nbhd_unique[match(neighborhood_names, nbhd_unique$neighborhood), , drop = FALSE]
X_nbhd <- cbind(
  `(Intercept)` = 1,
  walk_c = nbhd_unique$walk_score - mean(nbhd_unique$walk_score),
  transit_c = nbhd_unique$transit_score - mean(nbhd_unique$transit_score),
  bike_c = nbhd_unique$bike_score - mean(nbhd_unique$bike_score)
)
rownames(X_nbhd) <- neighborhood_names
p_nbhd <- ncol(X_nbhd)
nbhd_pred_names <- colnames(X_nbhd)
stopifnot(nrow(X_nbhd) == k)

message("airbnb neighborhood RE: n = ", l2, ", k = ", k, ", l1 = ", l1)
message("X_nbhd: ", nrow(X_nbhd), " x ", ncol(X_nbhd),
        " (", paste(nbhd_pred_names, collapse = ", "), ") — not passed to Block 1 yet")
message("Block 1: rindepNormalGamma_reg_multi (y = beta_mat ", k, " x ", l1, ")")
message("Block 2: rNormalGLM_reg_block_update")
message("Gibbs: n_burn = ", n_burn, ", n_sim = ", n_sim)

airbnb_dat$eta_proxy <- log(y_full + 1)
ps_glm <- Prior_Setup(
  eta_proxy ~ rating_c + room_type,
  family = gaussian(),
  data = airbnb_dat
)

hyper_Sigma_diag <- diag(diag(ps_glm$Sigma))
dimnames(hyper_Sigma_diag) <- list(beta_names, beta_names)
## Block 1 priors: fixed for entire Gibbs run (same prior_list_vector every iteration).
prior_list_vector <- build_prior_list_vector(ps_glm, beta_names)
prior_block2 <- list(
  mu = as.numeric(ps_glm$mu),
  Sigma = hyper_Sigma_diag,
  dispersion = 1,
  ddef = FALSE
)
fam <- poisson()

init_b2 <- rNormalGLM_reg_block_update(
  y = y_full, x = X_full, block = block_full,
  prior_list = prior_block2, family = fam,
  Gridtype = 2L, n_envopt = 1L,
  use_parallel = FALSE, use_opencl = use_opencl_random,
  verbose = FALSE, progbar = FALSE
)
beta_mat <- init_b2$coefficients
hyper_init <- block1_hyper_draw(
  beta_mat, prior_list_vector, beta_names,
  use_opencl_fixed = use_opencl_fixed
)
hyper_mu_loc <- hyper_init$mu
hyper_Sigma_loc <- hyper_init$Sigma

# --- Gibbs -------------------------------------------------------------------
run_gibbs <- function(store = FALSE) {
  hyper_mu_out <- if (store) matrix(0, n_sim, l1) else NULL
  hyper_disp_out <- if (store) matrix(0, n_sim, l1) else NULL
  coef_out <- if (store) {
    a <- array(0, c(n_sim, k, l1))
    dimnames(a) <- list(NULL, neighborhood_names, beta_names)
    a
  } else NULL

  beta_loc <- beta_mat
  hyper_mu_loc <- hyper_init$mu
  hyper_Sigma_loc <- hyper_init$Sigma
  t0 <- Sys.time()

  burn_time <- system.time({
    for (iter in seq_len(n_burn)) {
      b2 <- rNormalGLM_reg_block_update(
        y = y_full, x = X_full, block = block_full,
        prior_list = list(
          mu = hyper_mu_loc, Sigma = hyper_Sigma_loc,
          dispersion = 1, ddef = FALSE
        ),
        family = fam, Gridtype = 2L, n_envopt = 1L,
        use_parallel = FALSE, use_opencl = use_opencl_random,
        verbose = FALSE, progbar = FALSE
      )
      beta_loc <- b2$coefficients

      hyper_draw <- block1_hyper_draw(
        beta_loc, prior_list_vector, beta_names,
        use_opencl_fixed = use_opencl_fixed
      )
      hyper_mu_loc <- hyper_draw$mu
      hyper_Sigma_loc <- hyper_draw$Sigma

      if (iter == 1L) {
        sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        gibbs_message_first_iter_estimates("", sec, n_burn, n_sim)
      }
    }
  })

  report_every <- gibbs_report_interval(n_sim)
  t_main0 <- Sys.time()

  sim_time <- system.time({
    for (iter in seq_len(n_sim)) {
      b2 <- rNormalGLM_reg_block_update(
        y = y_full, x = X_full, block = block_full,
        prior_list = list(
          mu = hyper_mu_loc, Sigma = hyper_Sigma_loc,
          dispersion = 1, ddef = FALSE
        ),
        family = fam, Gridtype = 2L, n_envopt = 1L,
        use_parallel = FALSE, use_opencl = use_opencl_random,
        verbose = FALSE, progbar = FALSE
      )
      beta_loc <- b2$coefficients
      hyper_draw <- block1_hyper_draw(
        beta_loc, prior_list_vector, beta_names,
        use_opencl_fixed = use_opencl_fixed
      )
      hyper_mu_loc <- hyper_draw$mu
      hyper_Sigma_loc <- hyper_draw$Sigma
      if (store) {
        hyper_mu_out[iter, ] <- hyper_draw$mu
        hyper_disp_out[iter, ] <- hyper_draw$dispersion
        coef_out[iter, , ] <- beta_loc
      }
      if (iter %% report_every == 0L || iter == n_sim) {
        el <- difftime(Sys.time(), t_main0, units = "secs")
        message("Main ", iter, "/", n_sim, " — ", fmt_hms(el))
      }
    }
  })

  list(
    burn_time = burn_time,
    sim_time = sim_time,
    hyper_mu_out = hyper_mu_out,
    hyper_disp_out = hyper_disp_out,
    coef_block_out = coef_out,
    hyper_mu_final = hyper_mu_loc,
    hyper_Sigma_final = hyper_Sigma_loc,
    beta_final = beta_loc
  )
}

message("\n========== rindepNormalGamma_reg Block Gibbs ==========")
set.seed(123)
fit <- run_gibbs(store = TRUE)

t_total <- sum(as.numeric(fit$burn_time["elapsed"]), as.numeric(fit$sim_time["elapsed"]))
message("TOTAL: ", fmt_hms(t_total), " (", signif(t_total / n_gibbs, 4), " s/iter)")

benchmark <- list(
  hyper_sampler = "rindepNormalGamma_reg",
  block1 = "rindepNormalGamma_reg_multi",
  prior_list_vector = prior_list_vector,
  n = l2, k = k, l1 = l1,
  p_nbhd = p_nbhd,
  X_nbhd = X_nbhd,
  nbhd_pred_names = nbhd_pred_names,
  timing_seconds = list(total = t_total, per_iter = t_total / n_gibbs),
  hyper_mu_posterior_mean = colMeans(fit$hyper_mu_out),
  hyper_disp_posterior_mean = colMeans(fit$hyper_disp_out),
  beta_names = beta_names,
  timestamp = Sys.time()
)

coda_mu <- print_coda_hyper_chain(
  fit$hyper_mu_out, beta_names,
  "population-level mu (hyper mean across neighborhoods)"
)
if (!is.null(coda_mu)) {
  benchmark$coda_population_mu <- coda_mu
}

disp_names <- paste0("sigma2_", beta_names)
colnames(fit$hyper_disp_out) <- disp_names
message("\n--- Posterior mean dispersion (sigma^2) per random-effect dimension ---")
print(round(colMeans(fit$hyper_disp_out), 4))

coda_disp <- print_coda_hyper_chain(
  fit$hyper_disp_out, disp_names,
  "population-level dispersion sigma^2 (one per RE dimension, Block 1)"
)
if (!is.null(coda_disp)) {
  benchmark$coda_population_dispersion <- coda_disp
}

coef_mean <- apply(fit$coef_block_out, c(2L, 3L), mean)
avg_neighborhood <- colMeans(coef_mean)
if (!is.null(coda_mu)) {
  coda_means <- vapply(beta_names, function(nm) {
    coda_mu$summary[[1]][nm, "Mean"]
  }, numeric(1))
  cmp <- rbind(
    CODA_mu = coda_means,
    mean_neighborhood_beta = avg_neighborhood,
    diff = coda_means - avg_neighborhood
  )
  colnames(cmp) <- beta_names
  message("\n--- Population mu (CODA) vs mean of neighborhood posterior means ---")
  message("(Close values expected; not identical — mu is Block 1 draw, Average pools Block 2 betas)")
  print(round(cmp, 4))
  benchmark$population_mu_vs_neighborhood_average <- cmp
}
coef_table <- rbind(coef_mean, Average = colMeans(coef_mean))
rownames(coef_table) <- c(neighborhood_names, "Average")
colnames(coef_table) <- beta_names
message("\n--- Posterior mean beta by neighborhood ---")
print(round(coef_table, 4))
benchmark$coef_table <- coef_table

saveRDS(benchmark, file.path(root, "data-raw", OUT_RDS))
message("Wrote ", file.path(root, "data-raw", OUT_RDS))
message("Done.")

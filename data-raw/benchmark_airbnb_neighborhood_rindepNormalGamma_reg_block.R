# =============================================================================
# Airbnb neighborhood block Gibbs — Block 1: rindepNormalGamma_reg
# =============================================================================
# Hierarchical Poisson mixed model:
#   y_b | b_b ~ Poisson( Z_b  b_b )                        (likelihood)
#   b_b | beta, D ~ N( beta, D )                            (random effects prior)
#   beta_j, D_jj | b_mat ~ rindepNormalGamma_reg per j     (hyper / fixed effects)
#
# Notation follows lme4 / mixed-model convention:
#   Z    : listing-level design matrix  (n x p)  — random effects design
#   b    : neighborhood random effects  (k x p)
#   beta : population mean of b         (p)      — fixed effects
#   D    : diagonal covariance of b     (p x p)
#   X    : (not used here; see _covariates variant for X * beta hyper regression)
#
# Usage:
#   Rscript data-raw/benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block.R
#   Rscript data-raw/benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block.R quick

args      <- commandArgs(trailingOnly = TRUE)
run_quick <- any(tolower(args) == "quick")

pkgload::load_all(export_all = FALSE)
if (!requireNamespace("bayesrules", quietly = TRUE)) stop("Install bayesrules.")

OUT_RDS <- "data-raw/Airbnb_neighborhood_rindepNormalGamma_reg_benchmark.rds"

# --- Settings ----------------------------------------------------------------
n_burn <- if (run_quick) 5L   else 200L
n_sim  <- if (run_quick) 10L  else 1000L

# --- Data --------------------------------------------------------------------
data(list = if (run_quick) "airbnb_small" else "airbnb", package = "bayesrules")
airbnb_dat <- if (run_quick) airbnb_small else airbnb

airbnb_dat$rating_c  <- airbnb_dat$rating - mean(airbnb_dat$rating)
airbnb_dat$room_type <- factor(airbnb_dat$room_type)
airbnb_dat <- airbnb_dat[complete.cases(
  airbnb_dat[, c("reviews", "rating_c", "room_type", "neighborhood")]
), ]

Z          <- model.matrix(reviews ~ rating_c + room_type, data = airbnb_dat)
y          <- airbnb_dat$reviews
grp        <- factor(airbnb_dat$neighborhood)
n          <- nrow(Z)          # total observations
p          <- ncol(Z)          # number of random-effect dimensions (= ncol Z)
k          <- nlevels(grp)     # number of neighborhoods
z_names    <- colnames(Z)      # names of Z columns (= names of b dimensions)
grp_names  <- levels(grp)

message("n = ", n, "  k = ", k, "  p = ", p,
        "  n_burn = ", n_burn, "  n_sim = ", n_sim)

# --- Identifiability preflight -----------------------------------------------
id_check <- lmebayes:::block_check_identifiability_xy(Z, grp, on_failure = "stop")
stopifnot(id_check$action == "proceed")

# --- Priors ------------------------------------------------------------------
# Calibrate from a Gaussian proxy of the log-count surface
airbnb_dat$eta_proxy <- log(y + 1)
ps <- Prior_Setup(eta_proxy ~ rating_c + room_type, family = gaussian(), data = airbnb_dat)

# Block beta fixed prior: one dIndependent_Normal_Gamma per column of Z
prior_beta <- lapply(seq_len(p), function(j) {
  pl <- list(mu    = ps$mu[j],
             Sigma = matrix(diag(ps$Sigma)[j], 1, 1),
             shape = ps$shape_ING,
             rate  = ps$rate,
             max_disp_perc = 0.99)
  attr(pl, "Prior Type") <- "dIndependent_Normal_Gamma"
  pl
})

# Block b starting prior
prior_b <- list(mu = as.numeric(ps$mu), Sigma = diag(diag(ps$Sigma)),
                dispersion = 1, ddef = FALSE)

# --- Initialise --------------------------------------------------------------
# Initial draw of b (k x p matrix of neighborhood random effects)
b_mat <- rNormalGLM_reg_block_update(
  y = y, x = Z, block = grp,
  prior_list = prior_b, family = poisson(),
  Gridtype = 2L, n_envopt = 1L, use_parallel = FALSE
)$coefficients
colnames(b_mat) <- z_names

# Initial draw of beta and D from Block beta
beta_draw <- multi_rindepNormalGamma_reg(
  n = 1L, y = b_mat,
  x = matrix(1, k, 1, dimnames = list(NULL, "(Intercept)")),
  prior_list = prior_beta, family = gaussian(), use_parallel = FALSE
)
beta_loc <- vapply(beta_draw, function(f) f$coefficients[1, 1], numeric(1))
D_loc    <- diag(vapply(beta_draw, function(f) f$dispersion[1], numeric(1)))

# --- Gibbs loop --------------------------------------------------------------
beta_out <- matrix(0, n_sim, p, dimnames = list(NULL, z_names))
D_out    <- matrix(0, n_sim, p, dimnames = list(NULL, paste0("D_", z_names)))
b_out    <- array(0, c(n_sim, k, p), dimnames = list(NULL, grp_names, z_names))

set.seed(123)
t0 <- Sys.time()

one_iter <- function(beta, D) {
  # Block b: draw neighborhood random effects given current beta, D
  b_draw <- rNormalGLM_reg_block_update(
    y = y, x = Z, block = grp,
    prior_list = list(mu = beta, Sigma = D, dispersion = 1, ddef = FALSE),
    family = poisson(), Gridtype = 2L, n_envopt = 1L,
    use_parallel = FALSE, verbose = FALSE, progbar = FALSE
  )
  b_loc           <- b_draw$coefficients
  colnames(b_loc) <- z_names
  # Block beta: draw population mean and variance given current b
  beta_draw <- multi_rindepNormalGamma_reg(
    n = 1L, y = b_loc,
    x = matrix(1, k, 1, dimnames = list(NULL, "(Intercept)")),
    prior_list = prior_beta, family = gaussian(),
    use_parallel = FALSE, verbose = FALSE, progbar = FALSE
  )
  list(
    b    = b_loc,
    beta = vapply(beta_draw, function(f) f$coefficients[1, 1], numeric(1)),
    D    = diag(vapply(beta_draw, function(f) f$dispersion[1], numeric(1)))
  )
}

report_every_burn <- max(1L, n_burn %/% 10L)
for (iter in seq_len(n_burn)) {
  s <- one_iter(beta_loc, D_loc)
  beta_loc <- s$beta;  D_loc <- s$D
  if (iter %% report_every_burn == 0L || iter == n_burn)
    message("Burn-in ", iter, "/", n_burn,
            "  (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s)")
}

t1 <- Sys.time()
report_every_sim <- max(1L, n_sim %/% 10L)
for (iter in seq_len(n_sim)) {
  s <- one_iter(beta_loc, D_loc)
  beta_loc <- s$beta;  D_loc <- s$D
  beta_out[iter, ]   <- beta_loc
  D_out[iter, ]      <- diag(D_loc)
  b_out[iter, , ]    <- s$b
  if (iter %% report_every_sim == 0L || iter == n_sim)
    message("Sim ", iter, "/", n_sim,
            "  (", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), " s)")
}
t_sim <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
message("Simulation done (", round(t_sim, 1), " s  |  ", round(t_sim / n_sim, 4), " s/iter)")

# --- Results -----------------------------------------------------------------
message("\n--- Posterior means: beta (population mean of b) ---")
print(round(colMeans(beta_out), 4))

message("\n--- Posterior means: D (between-neighborhood variance) ---")
print(round(colMeans(D_out), 4))

b_mean    <- apply(b_out, c(2, 3), mean)
b_table   <- rbind(b_mean, Average = colMeans(b_mean))
message("\n--- Posterior mean b by neighborhood (and average) ---")
print(round(b_table, 4))

if (requireNamespace("coda", quietly = TRUE)) {
  as_mcmc <- function(x) {
    m <- coda::mcmc(as.matrix(x))
    if (length(dim(m)) == 3L &&
        exists("mcmcUpgrade", where = asNamespace("coda"), mode = "function")) {
      m <- coda::mcmcUpgrade(m)
    }
    m
  }
  message("\n--- CODA: beta ---")
  print(summary(as_mcmc(beta_out)))
  message("ESS:"); print(round(coda::effectiveSize(as_mcmc(beta_out)), 1))
}

benchmark <- list(
  hyper_sampler = "rindepNormalGamma_reg",
  n = n, k = k, p = p, z_names = z_names,
  timing_s = list(total = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  per_iter = t_sim / n_sim),
  beta_posterior_mean = colMeans(beta_out),
  D_posterior_mean    = colMeans(D_out),
  b_table             = b_table,
  timestamp           = Sys.time()
)
saveRDS(benchmark, OUT_RDS)
message("Wrote ", OUT_RDS, "\nDone.")

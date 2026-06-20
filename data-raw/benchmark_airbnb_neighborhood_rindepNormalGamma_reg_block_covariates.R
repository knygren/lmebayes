# =============================================================================
# Airbnb neighborhood block Gibbs — Block beta: rindepNormalGamma_reg + X
# =============================================================================
#   y_b | b_b ~ Poisson( Z_b  b_b )
#   b_b | beta, D ~ N( X_b %*% beta, D )     beta is (q x p); X_b is row b of X
#   beta, D | b ~ rindepNormalGamma_reg( b ~ X ) per column of Z
#
#   Z    : listing-level design     (n x p)
#   X    : neighborhood-level design (k x q)
#   b    : random effects           (k x p)
#   beta : fixed effects on X       (q x p)   — full Block-beta regression
#   X_beta : X %*% beta             (k x p)   — prior mean of b per neighborhood
#   D    : diag variance of b       (p x p)
#
# Usage:
#   Rscript data-raw/benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block_covariates.R
#   Rscript data-raw/benchmark_airbnb_neighborhood_rindepNormalGamma_reg_block_covariates.R quick

args      <- commandArgs(trailingOnly = TRUE)
run_quick <- any(tolower(args) == "quick")

pkgload::load_all(export_all = FALSE)
if (!requireNamespace("bayesrules", quietly = TRUE)) stop("Install bayesrules.")

OUT_RDS <- "data-raw/Airbnb_neighborhood_rindepNormalGamma_reg_benchmark_covariates.rds"

# --- Settings ----------------------------------------------------------------
n_burn <- if (run_quick) 5L   else 200L
n_sim  <- if (run_quick) 10L  else 1000L

# --- Data --------------------------------------------------------------------
data(list = if (run_quick) "airbnb_small" else "airbnb", package = "bayesrules")
airbnb_dat <- if (run_quick) airbnb_small else airbnb

airbnb_dat$rating_c  <- airbnb_dat$rating - mean(airbnb_dat$rating)
airbnb_dat$room_type <- factor(airbnb_dat$room_type)
airbnb_dat <- airbnb_dat[complete.cases(
  airbnb_dat[, c("reviews", "rating_c", "room_type", "neighborhood",
                 "walk_score", "transit_score", "bike_score")]
), ]

Z         <- model.matrix(reviews ~ rating_c + room_type, data = airbnb_dat)
y         <- airbnb_dat$reviews
grp       <- factor(airbnb_dat$neighborhood)
n         <- nrow(Z)
p         <- ncol(Z)
k         <- nlevels(grp)
z_names   <- colnames(Z)
grp_names <- levels(grp)

nbhd_rows <- airbnb_dat[!duplicated(airbnb_dat$neighborhood), ]
nbhd_rows <- nbhd_rows[match(grp_names, nbhd_rows$neighborhood), ]
X <- cbind(
  `(Intercept)` = 1,
  walk_c    = nbhd_rows$walk_score    - mean(nbhd_rows$walk_score),
  transit_c = nbhd_rows$transit_score - mean(nbhd_rows$transit_score),
  bike_c    = nbhd_rows$bike_score    - mean(nbhd_rows$bike_score)
)
rownames(X) <- grp_names
q       <- ncol(X)
x_names <- colnames(X)

message("n = ", n, "  k = ", k, "  p = ", p, "  q = ", q,
        "  n_burn = ", n_burn, "  n_sim = ", n_sim)

# --- Identifiability preflight -----------------------------------------------
id_check <- lmebayes:::block_check_identifiability_xy(Z, grp, X_nbhd = X, on_failure = "stop")
stopifnot(id_check$action == "proceed")

# --- Priors ------------------------------------------------------------------
airbnb_dat$eta_proxy <- log(y + 1)
ps <- Prior_Setup(eta_proxy ~ rating_c + room_type, family = gaussian(), data = airbnb_dat)
prior_b_init <- list(mu = as.numeric(ps$mu), Sigma = diag(diag(ps$Sigma)),
                     dispersion = 1, ddef = FALSE)

make_prior_block_beta <- function(b_mat) {
  X_df <- as.data.frame(X[, -1, drop = FALSE])
  lapply(seq_len(p), function(j) {
    X_df$b_j <- b_mat[, j]
    ps_j <- Prior_Setup(b_j ~ walk_c + transit_c + bike_c,
                        family = gaussian(), data = X_df)
    pl <- list(mu = as.numeric(ps_j$mu), Sigma = ps_j$Sigma,
               shape = ps_j$shape_ING, rate = ps_j$rate, max_disp_perc = 0.99)
    attr(pl, "Prior Type") <- "dIndependent_Normal_Gamma"
    pl
  })
}

# q x p coefficients from one Block-beta draw (rows = X, cols = Z)
beta_from_draw <- function(block_beta_draw) {
  b <- sapply(block_beta_draw, function(f) f$coefficients[1, ], simplify = "matrix")
  dimnames(b) <- list(x_names, z_names)
  b
}

# k x p prior means for b: X %*% beta
X_beta_from_beta <- function(beta) {
  Xb <- X %*% beta
  dimnames(Xb) <- list(grp_names, z_names)
  Xb
}

# --- Initialise --------------------------------------------------------------
b_mat <- rNormalGLM_reg_block_update(
  y = y, x = Z, block = grp,
  prior_list = prior_b_init, family = poisson(),
  Gridtype = 2L, n_envopt = 1L, use_parallel = FALSE
)$coefficients
colnames(b_mat) <- z_names

prior_block_beta <- make_prior_block_beta(b_mat)

block_beta_draw <- multi_rindepNormalGamma_reg(
  n = 1L, y = b_mat, x = X,
  prior_list = prior_block_beta, family = gaussian(), use_parallel = FALSE
)
beta_loc    <- beta_from_draw(block_beta_draw)
X_beta_loc  <- X_beta_from_beta(beta_loc)
D_loc       <- diag(vapply(block_beta_draw, function(f) f$dispersion[1], numeric(1)))

# --- Gibbs loop --------------------------------------------------------------
beta_out     <- array(0, c(n_sim, q, p), dimnames = list(NULL, x_names, z_names))
X_beta_out   <- array(0, c(n_sim, k, p), dimnames = list(NULL, grp_names, z_names))
D_out        <- matrix(0, n_sim, p, dimnames = list(NULL, paste0("D_", z_names)))
b_out        <- array(0, c(n_sim, k, p), dimnames = list(NULL, grp_names, z_names))

set.seed(123)
t0 <- Sys.time()

one_iter <- function(beta, D) {
  # Block b: b_b ~ N( (X %*% beta)_b, D )
  X_beta <- X_beta_from_beta(beta)
  b_draw <- rNormalGLM_reg_block_update(
    y = y, x = Z, block = grp,
    prior_list = list(mu = t(X_beta), Sigma = D, dispersion = 1, ddef = FALSE),
    family = poisson(), Gridtype = 2L, n_envopt = 1L,
    use_parallel = FALSE, verbose = FALSE, progbar = FALSE
  )
  b_loc           <- b_draw$coefficients
  colnames(b_loc) <- z_names
  # Block beta: b[,j] ~ X for each j
  block_beta_draw <- multi_rindepNormalGamma_reg(
    n = 1L, y = b_loc, x = X,
    prior_list = prior_block_beta, family = gaussian(),
    use_parallel = FALSE, verbose = FALSE, progbar = FALSE
  )
  beta_new <- beta_from_draw(block_beta_draw)
  list(
    b       = b_loc,
    beta    = beta_new,
    X_beta  = X_beta,
    D       = diag(vapply(block_beta_draw, function(f) f$dispersion[1], numeric(1)))
  )
}

report_every_burn <- max(1L, n_burn %/% 10L)
for (iter in seq_len(n_burn)) {
  s <- one_iter(beta_loc, D_loc)
  beta_loc <- s$beta
  X_beta_loc <- s$X_beta
  D_loc <- s$D
  if (iter %% report_every_burn == 0L || iter == n_burn)
    message("Burn-in ", iter, "/", n_burn,
            "  (", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s)")
}

t1 <- Sys.time()
report_every_sim <- max(1L, n_sim %/% 10L)
for (iter in seq_len(n_sim)) {
  s <- one_iter(beta_loc, D_loc)
  beta_loc <- s$beta
  X_beta_loc <- s$X_beta
  D_loc <- s$D
  beta_out[iter, , ]   <- beta_loc
  X_beta_out[iter, , ] <- s$X_beta
  D_out[iter, ]        <- diag(D_loc)
  b_out[iter, , ]       <- s$b
  if (iter %% report_every_sim == 0L || iter == n_sim)
    message("Sim ", iter, "/", n_sim,
            "  (", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), " s)")
}
t_sim <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
message("Simulation done (", round(t_sim, 1), " s  |  ", round(t_sim / n_sim, 4), " s/iter)")

# --- Results -----------------------------------------------------------------
beta_mean <- apply(beta_out, c(2, 3), mean)
dimnames(beta_mean) <- list(X = x_names, Z = z_names)

beta_sd <- apply(beta_out, c(2, 3), sd)
dimnames(beta_sd) <- list(X = x_names, Z = z_names)

message("\n--- beta posterior mean (q x p matrix: rows = X, cols = Z) ---")
print(round(beta_mean, 4))

message("\n--- beta posterior SD ---")
print(round(beta_sd, 4))

message("\n--- beta: walk / transit / bike rows only ---")
print(round(beta_mean[x_names != "(Intercept)", , drop = FALSE], 4))

message("\n--- Posterior means: X_beta = X %*% beta (prior mean of b, k x p) ---")
print(round(apply(X_beta_out, c(2, 3), mean), 4))

message("\n--- Posterior means: D ---")
print(round(colMeans(D_out), 4))

b_mean  <- apply(b_out, c(2, 3), mean)
b_table <- rbind(b_mean, Average = colMeans(b_mean))
message("\n--- Posterior mean b by neighborhood ---")
print(round(b_table, 4))

if (requireNamespace("coda", quietly = TRUE)) {
  # coda::mcmc() needs a 2D matrix; a 3D slice gives an obsolete mcmc object
  as_mcmc <- function(x) {
    m <- coda::mcmc(as.matrix(x))
    if (length(dim(m)) == 3L &&
        exists("mcmcUpgrade", where = asNamespace("coda"), mode = "function")) {
      m <- coda::mcmcUpgrade(m)
    }
    m
  }

  message("\n--- CODA: beta by column of Z (each block is q coefficients on X) ---")
  for (j in seq_len(p)) {
    message("\n  Z column: ", z_names[j])
    m <- as_mcmc(beta_out[, , j])
    colnames(m) <- x_names
    print(summary(m))
  }
  score_rows <- x_names[x_names != "(Intercept)"]
  message("\n--- ESS: walk / transit / bike (rows of beta; cols = Z) ---")
  ess_mat <- matrix(NA_real_, nrow = length(score_rows), ncol = p,
                    dimnames = list(score_rows, z_names))
  for (j in seq_len(p)) {
    m <- as_mcmc(beta_out[, , j])
    colnames(m) <- x_names
    ess_mat[, j] <- coda::effectiveSize(m)[score_rows]
  }
  print(round(ess_mat, 1))
}

benchmark <- list(
  hyper_sampler = "rindepNormalGamma_reg + X",
  n = n, k = k, p = p, q = q,
  z_names = z_names, x_names = x_names, X = X,
  timing_s = list(total = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  per_iter = t_sim / n_sim),
  beta_posterior_mean  = beta_mean,
  beta_posterior_sd    = beta_sd,
  X_beta_posterior_mean = apply(X_beta_out, c(2, 3), mean),
  D_posterior_mean     = colMeans(D_out),
  b_table              = b_table,
  timestamp            = Sys.time()
)
saveRDS(benchmark, OUT_RDS)
message("Wrote ", OUT_RDS, "\nDone.")

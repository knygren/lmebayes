# Benchmark rNormalGLM_reg_block on bayesrules::airbnb (Bayes Rules Ch. 18 style).
# Listing-level scalar random effects (x_one intercept per listing).
# For neighborhood-level multivariate coef vectors, see:
#   benchmark_airbnb_neighborhood_rNormalGLM_reg_block.R
#
# Default: Block 2 uses rNormalGLM_reg_block() (C++ block path). No per-observation
# R loops unless you append `legacy` for cross-checks afterward.
#
#   (1) Listing-level two-block Gibbs (population theta ~ rating + room_type).
#   (2) Grouped neighborhoods: one multivariate Block-2 update (l1 > 1, l2_b > 1).
#   (3) LEGACY only with `legacy`: rglmb loop (listing Gibbs) and/or
#       rNormal_reg loop per neighborhood (grouped check).
#
#   Rscript data-raw/benchmark_airbnb_rNormalGLM_reg_block.R
#   Rscript data-raw/benchmark_airbnb_rNormalGLM_reg_block.R quick
#   Rscript data-raw/benchmark_airbnb_rNormalGLM_reg_block.R legacy
#
# Default: n_burn = 200, n_sim = 1000. Append `quick` for smoke test only.
#
# Requires: pkgload, bayesrules (Suggests), optional coda for Gibbs summaries.

args <- commandArgs(trailingOnly = TRUE)
run_legacy <- any(tolower(args) %in% c("legacy", "--legacy", "-l"))
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))
path_args <- args[!tolower(args) %in% c(
  "legacy", "--legacy", "-l",
  "quick", "--quick", "-q"
)]

root <- if (length(path_args) >= 1L) {
  normalizePath(path_args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload, e.g. install.packages('pkgload')")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules, e.g. install.packages('bayesrules')")
}

pkgload::load_all(export_all = FALSE)

rglmb_update <- function(y, x, family = gaussian(), pfamily,
                         offset = NULL, weights = 1, Gridtype = 2L,
                         n_envopt = NULL, use_parallel = TRUE,
                         use_opencl = FALSE, verbose = FALSE,
                         mu_from_x = TRUE) {
  if (missing(pfamily)) stop("'pfamily' is required.", call. = FALSE)
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("nrow(x) must equal length(y).", call. = FALSE)
  out <- rglmb(n = 1L, y = y, x = x, family = family, pfamily = pfamily,
               offset = offset, weights = weights, Gridtype = as.integer(Gridtype),
               n_envopt = n_envopt, use_parallel = use_parallel,
               use_opencl = use_opencl, verbose = verbose)
  beta  <- as.vector(out$coefficients[1L, ])
  names(beta) <- colnames(out$coefficients)
  disp  <- out$dispersion[1L]
  mu_all <- if (isTRUE(mu_from_x)) as.vector(x %*% beta) else NULL
  list(beta = beta, dispersion = disp, sigma_theta_sq = disp,
       sigma_theta = sqrt(disp), mu_all = mu_all, rglmb = out)
}

summ_time <- function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
}

fmt_hms <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs) || secs < 0) {
    secs <- 0
  }
  h <- floor(secs / 3600)
  rem <- secs - h * 3600
  m <- floor(rem / 60)
  s <- rem - m * 60
  sprintf("%d h %d min %.2f s", h, m, s)
}

n_burn <- if (run_quick) 5L else 200L
n_sim  <- if (run_quick) 10L else 1000L
n_gibbs <- n_burn + n_sim
pct_listings_gibbs <- if (run_quick) 0.15 else 1

## --- Data (Bayes Rules: reviews ~ rating + room_type) ------------------------
data("airbnb", package = "bayesrules")

airbnb_dat <- airbnb
airbnb_dat$rating_c <- airbnb_dat$rating - mean(airbnb_dat$rating)
airbnb_dat$room_type <- factor(airbnb_dat$room_type)
airbnb_dat <- airbnb_dat[complete.cases(airbnb_dat[, c("reviews", "rating_c", "room_type", "neighborhood")]), ]

form_x <- reviews ~ rating_c + room_type
X_full <- model.matrix(form_x, data = airbnb_dat)
y_full <- airbnb_dat$reviews
block_full <- factor(airbnb_dat$neighborhood)
l2 <- length(y_full)
l1 <- ncol(X_full)
k <- nlevels(block_full)

block_info <- glmbayes:::normalize_block(block_full, l2)
stopifnot(block_info$k == k)

message("airbnb: n = ", l2, ", neighborhoods k = ", k,
        ", l1 = ", l1, " (", paste(colnames(X_full), collapse = ", "), ")")
message("listings per neighborhood: min = ", min(block_info$l2_blocks),
        ", median = ", median(block_info$l2_blocks),
        ", max = ", max(block_info$l2_blocks))
message("Block Gibbs: n_burn = ", n_burn, ", n_sim = ", n_sim,
        " (", n_gibbs, " full iterations); Block 2 = rNormalGLM_reg_block")
if (run_quick) {
  message("Quick mode: listing subset pct = ", pct_listings_gibbs)
}
if (!run_legacy) {
  message("Legacy R-loop cross-checks skipped (append 'legacy').")
}

airbnb_dat$eta_proxy <- log(y_full + 1)
ps_glm <- Prior_Setup(
  eta_proxy ~ rating_c + room_type,
  family = gaussian(),
  data = airbnb_dat
)
prior_template <- list(
  mu = as.numeric(ps_glm$mu),
  Sigma = ps_glm$Sigma,
  dispersion = 1,
  ddef = FALSE
)

fam <- poisson()

# =============================================================================
# (1) Listing-level two-block Gibbs — Block 2: rNormalGLM_reg_block
# =============================================================================

message("\n========== (1) Listing-level Gibbs (Block 2: rNormalGLM_reg_block) ==========")

set.seed(42)
n_gibbs_idx <- max(50L, round(pct_listings_gibbs * l2))
idx_gibbs <- if (pct_listings_gibbs < 1) {
  sample(l2, size = n_gibbs_idx)
} else {
  seq_len(l2)
}
n_gibbs_idx <- length(idx_gibbs)

Bike_gibbs <- airbnb_dat[idx_gibbs, , drop = FALSE]
y_train <- y_full[idx_gibbs]
X_train <- X_full[idx_gibbs, , drop = FALSE]
n_train <- length(y_train)
p <- ncol(X_train)

theta <- log(y_train + 1)
data_pop <- data.frame(theta = theta, Bike_gibbs)
form_pop <- theta ~ rating_c + room_type
ps_pop <- Prior_Setup(form_pop, family = gaussian(), data = data_pop)
pfamily_pop <- dNormal_Gamma(
  ps_pop$mu, Sigma_0 = ps_pop$Sigma_0,
  ps_pop$shape, ps_pop$rate
)

x_one <- matrix(1, n_train, 1)
colnames(x_one) <- "(Intercept)"

block2_theta_rglmb_loop <- function(mu_all, sigma_theta_sq, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  theta_new <- numeric(n_train)
  for (i in seq_len(n_train)) {
    theta_new[i] <- rglmb(
      1L, y = y_train[i], x = matrix(1, 1, 1),
      family = fam,
      pfamily = dNormal(mu = mu_all[i], Sigma = sigma_theta_sq),
      Gridtype = 2L,
      use_parallel = FALSE,
      use_opencl = FALSE,
      verbose = FALSE
    )$coefficients[1, 1]
  }
  theta_new
}

gibbs_report_interval <- function(n) {
  max(1L, min(50L, as.integer(n %/% 10L)))
}

gibbs_message_first_iter_estimates <- function(tag, sec_per_iter, n_burn, n_sim) {
  message(tag, "After iteration 1 — estimated burn-in: ",
          fmt_hms(sec_per_iter * n_burn),
          "; estimated total simulation (burn-in + main): ",
          fmt_hms(sec_per_iter * (n_burn + n_sim)),
          " (", signif(sec_per_iter, 4), " s/iteration)")
}

run_block_gibbs <- function(store = FALSE, label = "") {
  beta_out <- if (store) matrix(0, nrow = n_sim, ncol = p) else NULL
  sigma_out <- if (store) numeric(n_sim) else NULL
  theta_loc <- theta
  tag <- if (nzchar(label)) paste0("[", label, "] ") else ""
  t_chain_start <- Sys.time()

  burn_time <- system.time({
    for (iter in seq_len(n_burn)) {
      b1 <- rglmb_update(
        y = theta_loc,
        x = X_train,
        family = gaussian(),
        pfamily = pfamily_pop,
        use_parallel = FALSE,
        verbose = FALSE
      )
      theta_loc <- rNormalGLM_reg_block_update(
        mu_all = b1$mu_all,
        sigma_theta_sq = b1$sigma_theta_sq,
        y = y_train,
        x = x_one,
        block = seq_len(n_train),
        family = fam,
        Gridtype = 2L,
        n_envopt = 1L,
        use_parallel = FALSE,
        use_opencl = FALSE,
        verbose = FALSE,
        progbar = FALSE
      )$theta
      if (iter == 1L) {
        sec_per_iter <- as.numeric(difftime(Sys.time(), t_chain_start, units = "secs"))
        gibbs_message_first_iter_estimates(tag, sec_per_iter, n_burn, n_sim)
      }
    }
  })
  t_burn_elapsed <- as.numeric(burn_time["elapsed"])
  sec_per_iter <- t_burn_elapsed / n_burn
  message(tag, "Burn-in complete: ", fmt_hms(t_burn_elapsed),
          " (", n_burn, " iterations, ",
          signif(sec_per_iter, 4), " s/iteration)")
  message(tag, "Estimated time remaining (main phase, ",
          n_sim, " iterations): ",
          fmt_hms(sec_per_iter * n_sim))

  report_every <- gibbs_report_interval(n_sim)
  t_main_start <- Sys.time()
  sim_time <- system.time({
    for (iter in seq_len(n_sim)) {
      b1 <- rglmb_update(
        y = theta_loc,
        x = X_train,
        family = gaussian(),
        pfamily = pfamily_pop,
        use_parallel = FALSE,
        verbose = FALSE
      )
      theta_loc <- rNormalGLM_reg_block_update(
        mu_all = b1$mu_all,
        sigma_theta_sq = b1$sigma_theta_sq,
        y = y_train,
        x = x_one,
        block = seq_len(n_train),
        family = fam,
        Gridtype = 2L,
        n_envopt = 1L,
        use_parallel = FALSE,
        use_opencl = FALSE,
        verbose = FALSE,
        progbar = FALSE
      )$theta
      if (store) {
        beta_out[iter, ] <- b1$beta
        sigma_out[iter] <- b1$sigma_theta
      }
      if (iter == 1L) {
        sec_main <- as.numeric(difftime(Sys.time(), t_main_start, units = "secs"))
        message(tag, "After 1 main iteration — estimated main phase: ",
                fmt_hms(sec_main * n_sim),
                "; estimated total simulation: ",
                fmt_hms(t_burn_elapsed + sec_main * n_sim),
                " (", signif(sec_main, 4), " s/iteration in main)")
      }
      if (iter %% report_every == 0L || iter == n_sim) {
        elapsed_main <- as.numeric(difftime(Sys.time(), t_main_start, units = "secs"))
        remaining <- if (iter > 0L) {
          elapsed_main / iter * (n_sim - iter)
        } else {
          sec_per_iter * n_sim
        }
        message(tag, "Main: ", iter, "/", n_sim,
                " — elapsed ", fmt_hms(elapsed_main),
                ", ETA ", fmt_hms(remaining))
      }
    }
  })

  list(
    burn_time = burn_time,
    sim_time = sim_time,
    beta_out = beta_out,
    sigma_out = sigma_out,
    theta_final = theta_loc
  )
}

message("Listing subset n = ", n_train, ", p = ", p,
        "; Block 2 blocks k = n, l1 = 1 per listing")
message("Started: ", format(Sys.time(), usetz = TRUE))

set.seed(123)
gibbs_blk <- run_block_gibbs(store = TRUE)

t_burn <- as.numeric(gibbs_blk$burn_time["elapsed"])
t_sim  <- as.numeric(gibbs_blk$sim_time["elapsed"])
t_total <- t_burn + t_sim

message("\n--- timing: listing-level Block Gibbs ---")
message("burn-in (", n_burn, "): ", fmt_hms(t_burn), " (", signif(t_burn, 4), " s)")
message("main (", n_sim, "):     ", fmt_hms(t_sim), " (", signif(t_sim, 4), " s)")
message("TOTAL (", n_gibbs, "):   ", fmt_hms(t_total), " (", signif(t_total, 4), " s)")
message("mean seconds per iteration: ", signif(t_total / n_gibbs, 4))

benchmark <- list(
  dataset = "bayesrules::airbnb",
  outcome = "reviews",
  n_full = l2,
  k_neighborhoods = k,
  l1_grouped = l1,
  quick_mode = run_quick,
  block2_method = "rNormalGLM_reg_block",
  listing_gibbs = list(
    n_train = n_train,
    p = p,
    n_burn = n_burn,
    n_sim = n_sim,
    timing_seconds = list(
      burn_in = t_burn,
      main = t_sim,
      total = t_total,
      per_iteration_mean = t_total / n_gibbs
    ),
    population_formula = deparse(form_pop, width.cutoff = 500L),
    posterior_mean_beta = colMeans(gibbs_blk$beta_out),
    posterior_mean_sigma_theta = mean(gibbs_blk$sigma_out)
  ),
  timestamp = Sys.time()
)

if (requireNamespace("coda", quietly = TRUE) && n_sim >= 20L) {
  beta_names <- colnames(X_train)
  mcmc_gibbs <- coda::mcmc(cbind(gibbs_blk$beta_out, sigma_theta = gibbs_blk$sigma_out))
  colnames(mcmc_gibbs) <- c(beta_names, "sigma_theta")
  message("\n--- CODA summary (listing-level population block) ---")
  print(summary(mcmc_gibbs))
  benchmark$listing_gibbs$effective_size <- as.list(coda::effectiveSize(mcmc_gibbs))
}

# =============================================================================
# (2) Grouped neighborhoods — one Block-2 update (multivariate blocks)
# =============================================================================

message("\n========== (2) Grouped blocks: one rNormalGLM_reg_block update ==========")
message("Started: ", format(Sys.time(), usetz = TRUE))

set.seed(2026)
n_time_grp <- if (run_quick) 2L else 3L
time_blk_grp <- numeric(n_time_grp)
for (t in seq_len(n_time_grp)) {
  time_blk_grp[t] <- system.time({
    out_grp <- rNormalGLM_reg_block(
      n = 1L,
      y = y_full,
      x = X_full,
      block = block_full,
      prior_list = prior_template,
      family = fam,
      Gridtype = 2L,
      n_envopt = 1L,
      use_parallel = FALSE,
      use_opencl = FALSE,
      verbose = FALSE,
      progbar = FALSE
    )
  })["elapsed"]
}
s_blk_grp <- summ_time(time_blk_grp)
message("rNormalGLM_reg_block (k = ", k, ", l1 = ", l1, "):")
print(round(s_blk_grp, 3))

benchmark$grouped_one_update <- list(
  timing_seconds = list(rNormalGLM_reg_block = s_blk_grp),
  coef_mode = out_grp$coef.mode
)

# =============================================================================
# (3) LEGACY: R-loop cross-checks only (append `legacy` on command line)
# =============================================================================

if (run_legacy) {
  message("\n========== (3) LEGACY cross-checks (R loops) ==========")

  message("\n--- (3a) Listing Gibbs: Block 2 via per-listing rglmb loop ---")
  message("Started: ", format(Sys.time(), usetz = TRUE))
  set.seed(123)
  run_block_gibbs_legacy <- function(store = FALSE, label = "legacy") {
    beta_out <- if (store) matrix(0, nrow = n_sim, ncol = p) else NULL
    sigma_out <- if (store) numeric(n_sim) else NULL
    theta_loc <- theta
    tag <- if (nzchar(label)) paste0("[", label, "] ") else ""
    burn_time <- system.time({
      for (iter in seq_len(n_burn)) {
        b1 <- rglmb_update(
          y = theta_loc, x = X_train, family = gaussian(),
          pfamily = pfamily_pop, use_parallel = FALSE, verbose = FALSE
        )
        theta_loc <- block2_theta_rglmb_loop(b1$mu_all, b1$sigma_theta_sq)
      }
    })
    sim_time <- system.time({
      for (iter in seq_len(n_sim)) {
        b1 <- rglmb_update(
          y = theta_loc, x = X_train, family = gaussian(),
          pfamily = pfamily_pop, use_parallel = FALSE, verbose = FALSE
        )
        theta_loc <- block2_theta_rglmb_loop(b1$mu_all, b1$sigma_theta_sq)
        if (store) {
          beta_out[iter, ] <- b1$beta
          sigma_out[iter] <- b1$sigma_theta
        }
      }
    })
    list(
      burn_time = burn_time, sim_time = sim_time,
      beta_out = beta_out, sigma_out = sigma_out,
      theta_final = theta_loc
    )
  }

  gibbs_rglmb <- run_block_gibbs_legacy(store = FALSE)

  t_burn_l <- as.numeric(gibbs_rglmb$burn_time["elapsed"])
  t_sim_l  <- as.numeric(gibbs_rglmb$sim_time["elapsed"])
  t_total_l <- t_burn_l + t_sim_l

  message("burn-in: ", fmt_hms(t_burn_l), " (", signif(t_burn_l, 4), " s)")
  message("main:    ", fmt_hms(t_sim_l), " (", signif(t_sim_l, 4), " s)")
  message("TOTAL:   ", fmt_hms(t_total_l), " (", signif(t_total_l, 4), " s)")
  message("speedup block vs legacy (total): ",
          signif(t_total_l / t_total, 3), "x")

  message("\n--- one Block-2 update (listing level), for reference ---")
  b1_ref <- rglmb_update(
    y = gibbs_blk$theta_final,
    x = X_train,
    family = gaussian(),
    pfamily = pfamily_pop,
    use_parallel = FALSE,
    verbose = FALSE
  )
  n_time_b2 <- if (run_quick) 2L else 3L
  time_blk2 <- numeric(n_time_b2)
  time_rglmb2 <- numeric(n_time_b2)
  for (t in seq_len(n_time_b2)) {
    time_blk2[t] <- system.time({
      rNormalGLM_reg_block_update(
        mu_all = b1_ref$mu_all,
        sigma_theta_sq = b1_ref$sigma_theta_sq,
        y = y_train,
        x = x_one,
        block = seq_len(n_train),
        family = fam,
        Gridtype = 2L,
        n_envopt = 1L,
        use_parallel = FALSE,
        use_opencl = FALSE,
        verbose = FALSE,
        progbar = FALSE
      )
    })["elapsed"]
    time_rglmb2[t] <- system.time({
      block2_theta_rglmb_loop(b1_ref$mu_all, b1_ref$sigma_theta_sq)
    })["elapsed"]
  }
  s_blk2 <- summ_time(time_blk2)
  s_rglmb2 <- summ_time(time_rglmb2)
  message("rNormalGLM_reg_block (one Block-2 step):")
  print(round(s_blk2, 3))
  message("rglmb loop (one Block-2 step):")
  print(round(s_rglmb2, 3))

  benchmark$legacy_listing_gibbs <- list(
    block2_method = "rglmb_loop",
    timing_seconds = list(
      burn_in = t_burn_l,
      main = t_sim_l,
      total = t_total_l,
      per_iteration_mean = t_total_l / n_gibbs,
      block2_only_rNormalGLM_reg_block = s_blk2,
      block2_only_rglmb_loop = s_rglmb2
    ),
    speedup_total_block_vs_legacy = t_total_l / t_total
  )

  message("\n--- (3b) Grouped: rNormal_reg loop over neighborhoods ---")
  prior_block <- glmbayes:::normalize_prior_for_blocks(
    prior_list = prior_template,
    prior_lists = NULL,
    block_info = block_info,
    l1 = l1
  )
  time_rnr_grp <- numeric(n_time_grp)
  for (t in seq_len(n_time_grp)) {
    time_rnr_grp[t] <- system.time({
      coef_rnr <- matrix(NA_real_, nrow = k, ncol = l1)
      for (b in seq_len(k)) {
        rows_b <- block_info$rows[[b]]
        out_b <- rNormal_reg(
          n = 1L,
          y = y_full[rows_b],
          x = X_full[rows_b, , drop = FALSE],
          prior_list = prior_block[[b]],
          family = fam,
          Gridtype = 2L,
          n_envopt = 1L,
          use_parallel = FALSE,
          use_opencl = FALSE,
          verbose = FALSE,
          progbar = FALSE
        )
        cb <- out_b$coefficients
        coef_rnr[b, ] <- if (is.matrix(cb)) cb[1L, ] else as.numeric(cb)
      }
    })["elapsed"]
  }
  s_rnr_grp <- summ_time(time_rnr_grp)
  message("rNormal_reg per neighborhood:")
  print(round(s_rnr_grp, 3))
  message("speedup grouped block vs legacy loop (mean): ",
          signif(s_rnr_grp["mean"] / s_blk_grp["mean"], 3), "x")
  max_mode_diff <- max(abs(coef_rnr - out_grp$coef.mode))
  message("max |mode_rNormal_reg_loop - mode_block| = ", signif(max_mode_diff, 6))

  benchmark$legacy_grouped <- list(
    timing_seconds = list(rNormal_reg_loop = s_rnr_grp),
    speedup_block_vs_loop = s_rnr_grp["mean"] / s_blk_grp["mean"],
    max_mode_abs_diff = max_mode_diff
  )
}

out_path <- file.path(root, "data-raw", "Airbnb_block_reg_benchmark.rds")
saveRDS(benchmark, out_path)
message("\nWrote summary: ", out_path)
message("Done.")

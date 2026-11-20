## Ridge weight, part 2, Poisson: groups draws under the full ridge model.
##
## The ridge model replaces each group log-likelihood l_j by
##   (1 - rwt) l_j + rwt T_j,
## where T_j is the second-order Taylor polynomial of l_j at the reference
## point m_j (the ICM mode). The sweep scales the weights by 1 - rwt, and
## .two_block_ridge_attach() puts rwt T_j into each group's prior. This
## script checks, with the population coefficients fixed at the ICM start:
##   1. the ridge curvature C_j (central differences of f3) against the
##      exact Poisson curvature X' diag(lambda) X at m_j;
##   2. the compiled groups draws against 2-D quadrature of
##        prior(b_j) * exp((1 - rwt) l_j(b_j) + rwt T_j(b_j));
##   3. the R and C++ sweeps with every draw replaced by its posterior mode,
##      ridge attached (these must agree exactly).
##
## Model: airbnb_small, reviews ~ walk_c + rating_c + (1 + rating_c ||
## neighborhood), as in check_ridge_part1_poisson.R.
## Draws are iid (a fresh envelope draw per call), so the SE of a mean is
## sqrt(var / n_draws) and the SE of a variance ratio is about
## sqrt(2 / (n_draws - 1)).
##
##   Rscript C:/Rpackages/lmebayes/data-raw/check_ridge_part2_poisson.R

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("C:/Rpackages/lmebayesCore", quiet = TRUE)
ns <- asNamespace("lmebayesCore")
get_ns <- function(name) get(name, envir = ns)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

n_draws <- 2000L
rwt_values <- c(0.25, 0.75)

data(airbnb_small, package = "bayesrules")
dat <- airbnb_small
dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
dat$walk_c <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
dat <- dat[stats::complete.cases(dat[, c(
  "reviews", "rating_c", "walk_c", "neighborhood"
)]), ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))
form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

design <- lmebayesCore::model_setup(form, data = dat, family = poisson())
ps <- lmebayesCore::Prior_Setup_GLMM(form, data = dat, family = poisson(),
                                     pop.pwt = 0.01)
pf <- lmebayesCore::pfamily_list(ps)
re <- design$groupef.names
gl <- levels(design$group)
block1_prior <- list(ddef = TRUE, P = get_ns(".rLMM_P_from_pfamily_list")(pf, re))

icm <- get_ns(".two_block_icm_at_start")(
  design = design, prior_list = block1_prior, pfamily_list = pf,
  re_names = re, family = poisson(), tol = 1e-10, maxit = 200L
)
map <- get_ns(".two_block_index_map")(design, re, gl, icm$start)
spec <- get_ns(".two_block_spec")(design, block1_prior, pf, poisson(), re, map)
state <- get_ns(".two_block_init")(1L, icm$start, icm$b_start, spec, design,
                                   re, gl)
mu_D <- get_ns(".two_block_groups_prior_mean")(spec, state$fixef, 1L)
attach_ridge <- get_ns(".two_block_ridge_attach")

D <- spec$D
l1 <- ncol(D)
b_ref <- as.matrix(icm$b_start)
if (!is.null(rownames(b_ref))) b_ref <- b_ref[gl, , drop = FALSE]
if (!is.null(colnames(b_ref))) b_ref <- b_ref[, re, drop = FALSE]
ref_point <- function(j) {
  m <- numeric(l1)
  m[spec$d_col] <- b_ref[j, ]
  m
}

offset_r <- function(r) {
  if (length(spec$offset1) == 0L) 0 else spec$offset1[r]
}

## Exact Poisson pieces at the reference point: log-likelihood, score,
## curvature.
pois_terms <- function(j, m) {
  r <- spec$obs_rows[[j]]
  Xr <- D[r, , drop = FALSE]
  eta <- as.numeric(Xr %*% m) + offset_r(r)
  lam <- exp(eta)
  w <- spec$wt1[r]
  list(
    l = sum(w * (spec$y[r] * eta - lam)),
    g = as.numeric(crossprod(Xr, w * (spec$y[r] - lam))),
    C = crossprod(Xr, Xr * (w * lam))
  )
}

## --- 1. curvature cross-check ----------------------------------------------
cat("1. Ridge curvature vs exact Poisson curvature at the reference point\n")
spec_chk <- attach_ridge(spec, b_ref, 0.25, re, gl)
gap_C <- 0
for (j in seq_along(gl)) {
  C_ridge <- (spec_chk$ridge_P[[j]] - spec$P_D) / 0.25
  C_exact <- pois_terms(j, ref_point(j))$C
  gap_C <- max(gap_C, max(abs(C_ridge - C_exact)) / max(1, max(abs(C_exact))))
}
cat(sprintf("   max relative |C_ridge - C_exact| over %d groups = %.3g\n",
            length(gl), gap_C))

## --- 2. draws vs quadrature --------------------------------------------------
rows_n <- vapply(spec$obs_rows, length, integer(1))
pick <- c(
  rich   = which.max(rows_n),
  median = order(rows_n)[ceiling(length(rows_n) / 2)],
  sparse = which.min(rows_n)
)
pick <- pick[!duplicated(pick)]
cat("\n2. Neighborhoods checked:\n")
for (nm in names(pick)) {
  j <- pick[[nm]]
  r <- spec$obs_rows[[j]]
  cat(sprintf("  %-7s %-28s listings = %3d, reviews = %d\n",
              nm, gl[j], rows_n[j], as.integer(sum(spec$y[r]))))
}

log_target <- function(b, j, w, m, tm) {
  r <- spec$obs_rows[[j]]
  eta <- as.numeric(D[r, , drop = FALSE] %*% b) + offset_r(r)
  ll <- sum(spec$wt1[r] * (spec$y[r] * eta - exp(eta)))
  dm <- b - m
  taylor <- tm$l + sum(tm$g * dm) - 0.5 * sum(dm * (tm$C %*% dm))
  d <- b - mu_D[, j]
  (1 - w) * ll + w * taylor - 0.5 * sum(d * (spec$P_D %*% d))
}

## Quadrature on a grid of +/- 8 Laplace SDs around the mode.
quad_moments <- function(j, w, mode, n_grid = 241L) {
  m <- ref_point(j)
  tm <- pois_terms(j, m)
  r <- spec$obs_rows[[j]]
  Xr <- D[r, , drop = FALSE]
  lam <- exp(as.numeric(Xr %*% mode) + offset_r(r))
  H <- (1 - w) * crossprod(Xr, Xr * (spec$wt1[r] * lam)) + w * tm$C + spec$P_D
  s <- sqrt(diag(solve(H)))
  g1 <- seq(mode[1] - 8 * s[1], mode[1] + 8 * s[1], length.out = n_grid)
  g2 <- seq(mode[2] - 8 * s[2], mode[2] + 8 * s[2], length.out = n_grid)
  lp <- outer(seq_along(g1), seq_along(g2), Vectorize(function(a, b) {
    log_target(c(g1[a], g2[b]), j, w, m, tm)
  }))
  wq <- exp(lp - max(lp))
  wq <- wq / sum(wq)
  m1 <- sum(wq * g1[row(wq)])
  m2 <- sum(wq * g2[col(wq)])
  v1 <- sum(wq * (g1[row(wq)] - m1)^2)
  v2 <- sum(wq * (g2[col(wq)] - m2)^2)
  list(mean = c(m1, m2), var = c(v1, v2))
}

draw_cpp <- get_ns(".two_block_groups_draw_cpp")
cn <- colnames(design$D)
rows <- list()
for (w_rwt in rwt_values) {
  spec_w <- attach_ridge(spec, b_ref, w_rwt, re, gl)
  spec_w$wt1 <- spec_w$wt1 * (1 - w_rwt)
  t0 <- proc.time()[["elapsed"]]
  draws <- array(NA_real_, c(n_draws, ncol(D), length(pick)))
  modes <- NULL
  for (s in seq_len(n_draws)) {
    d <- draw_cpp(spec_w, mu_D)
    draws[s, , ] <- t(d$coef[pick, , drop = FALSE])
    if (is.null(modes)) modes <- d$coef_mode[pick, , drop = FALSE]
  }
  cat(sprintf("\nrwt = %g: %d draws in %.1f s\n", w_rwt, n_draws,
              proc.time()[["elapsed"]] - t0))
  for (a in seq_along(pick)) {
    j <- pick[[a]]
    q <- quad_moments(j, w_rwt, as.numeric(modes[a, ]))
    for (c in seq_len(ncol(D))) {
      x <- draws[, c, a]
      z <- (mean(x) - q$mean[c]) / sqrt(q$var[c] / n_draws)
      ratio <- stats::var(x) / q$var[c]
      z_var <- (ratio - 1) / sqrt(2 / (n_draws - 1))
      rows[[length(rows) + 1L]] <- data.frame(
        rwt = w_rwt, nbhd = gl[j], which = names(pick)[a], coef = cn[c],
        quad_mean = q$mean[c], draw_mean = mean(x), z_mean = z,
        quad_sd = sqrt(q$var[c]), var_ratio = ratio, z_var = z_var,
        stringsAsFactors = FALSE
      )
    }
  }
}
res <- do.call(rbind, rows)
cat("\nDraws vs quadrature (z for means; var_ratio = draw var / quad var)\n")
print(res, digits = 3, row.names = FALSE)

z_all <- c(res$z_mean, res$z_var)
cat(sprintf(
  "\n%d z-scores: max |z| = %.2f, |z| > 3: %d (about %.2f expected by chance)\n",
  length(z_all), max(abs(z_all)), sum(abs(z_all) > 3),
  length(z_all) * 2 * stats::pnorm(-3)
))

## --- 3. R vs C++ sweep, ridge attached, use_mode ---------------------------
n_chains <- 2L
n_sweeps <- 3L
start2 <- lapply(icm$start, function(g) g)
state2 <- get_ns(".two_block_init")(n_chains, start2, icm$b_start, spec,
                                    design, re, gl)
cat(sprintf("\n3. R vs C++ sweep at rwt = 0.25, %d chains, %d sweeps, use_mode\n",
            n_chains, n_sweeps))
spec_ridge <- attach_ridge(spec, b_ref, 0.25, re, gl)
spec_ridge_r <- spec_ridge
spec_ridge_r$wt1 <- spec_ridge$wt1 * (1 - 0.25)
sw_r <- get_ns(".two_block_sweep_r")(spec_ridge_r, state2, n_sweeps,
                                     use_mode = TRUE)
sw_c <- get_ns(".two_block_sweep_cpp")(
  spec_ridge, state2$b, state2$fixef, state2$iters, state2$iters_ranef,
  n_sweeps, TRUE, FALSE, "", 0.25
)
gap <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
gaps <- c(
  b = gap(sw_r$b, sw_c$b),
  fixef = max(vapply(seq_along(re), function(k) {
    gap(sw_r$fixef[[k]], sw_c$fixef[[k]])
  }, numeric(1))),
  sweep_mean = gap(sw_r$sweep_mean, sw_c$sweep_mean)
)
for (nm in names(gaps)) {
  cat(sprintf("   %-11s max |R - C++| = %.3g %s\n", nm, gaps[[nm]],
              if (gaps[[nm]] == 0) "exact" else "DIFFERS"))
}
sw_plain <- get_ns(".two_block_sweep_cpp")(
  spec, state2$b, state2$fixef, state2$iters, state2$iters_ranef,
  n_sweeps, TRUE, FALSE, "", 0.25
)
cat(sprintf("   ridge moves the group modes: max |b_ridge - b_part1| = %.3g\n",
            gap(sw_c$b, sw_plain$b)))

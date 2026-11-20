## Ridge weight, part 1: groups draws under likelihood^(1 - rwt).
##
## The C++ sweep draws the groups block with observation weights
## wt1 * (1 - rwt) (check_port_steps.R confirms the sweep passes exactly
## those weights). This script checks the other half: given those weights,
## the compiled groups sampler draws from
##   prior(b_j) * L_j(b_j)^(1 - rwt)
## for each state j. The reference is 2-D quadrature of that density.
##
## Model: demo Ex_17 (binomial logit, (1 + violent_i || state)), as in
## check_port_steps.R. The population coefficients are fixed at the ICM
## start, so each state's conditional prior is fixed.
## rwt = 0.25 and 0.75 (the target needs rwt > 0).
## Draws are iid (a fresh envelope draw per call), so the SE of a mean is
## sqrt(var / n_draws) and the SE of a variance ratio is about
## sqrt(2 / (n_draws - 1)).
##
##   Rscript C:/Rpackages/lmebayes/data-raw/check_ridge_part1.R

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("C:/Rpackages/lmebayesCore", quiet = TRUE)
ns <- asNamespace("lmebayesCore")
get_ns <- function(name) get(name, envir = ns)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

n_draws <- 2000L
rwt_values <- c(0.25, 0.75)

data(book_banning, package = "bayesrules")
dat <- book_banning[, c("state", "removed", "violent", "median_income")]
dat <- dat[stats::complete.cases(dat), ]
dat$removed_i <- as.integer(dat$removed == 1L | dat$removed == "1")
dat$violent_i <- as.integer(
  dat$violent == TRUE | dat$violent == 1L | dat$violent == "TRUE"
)
dat$income_c <- as.numeric(scale(dat$median_income))
form <- removed_i ~ violent_i + income_c + income_c:violent_i +
  (1 + violent_i || state)

design <- lmebayesCore::model_setup(form, data = dat, family = binomial())
ps <- lmebayesCore::Prior_Setup_GLMM(form, data = dat, family = binomial(),
                                     pop.pwt = 0.01)
pf <- lmebayesCore::pfamily_list(ps)
re <- design$groupef.names
gl <- levels(design$group)
block1_prior <- list(ddef = TRUE, P = get_ns(".rLMM_P_from_pfamily_list")(pf, re))

icm <- get_ns(".two_block_icm_at_start")(
  design = design, prior_list = block1_prior, pfamily_list = pf,
  re_names = re, family = binomial(), tol = 1e-10, maxit = 200L
)
map <- get_ns(".two_block_index_map")(design, re, gl, icm$start)
spec <- get_ns(".two_block_spec")(design, block1_prior, pf, binomial(), re, map)
state <- get_ns(".two_block_init")(1L, icm$start, icm$b_start, spec, design,
                                   re, gl)
mu_D <- get_ns(".two_block_groups_prior_mean")(spec, state$fixef, 1L)

## Three states: most observations, a state with no violent books (flat
## likelihood in the slope direction) when there is one, and a median one.
D <- spec$D
viol_col <- spec$d_col[match("violent_i", re)]
rows_n <- vapply(spec$obs_rows, length, integer(1))
n_viol <- vapply(spec$obs_rows, function(r) sum(D[r, viol_col]), numeric(1))
pick <- c(rich = which.max(rows_n))
if (any(n_viol == 0)) {
  pick <- c(pick, no_violent = which(n_viol == 0)[which.max(rows_n[n_viol == 0])])
}
pick <- c(pick, median = order(rows_n)[ceiling(length(rows_n) / 2)])
pick <- pick[!duplicated(pick)]
cat("States checked:\n")
for (nm in names(pick)) {
  j <- pick[[nm]]
  cat(sprintf("  %-11s %-4s rows = %3d, violent books = %d\n",
              nm, gl[j], rows_n[j], as.integer(n_viol[j])))
}

log_target <- function(b, j, w) {
  r <- spec$obs_rows[[j]]
  eta <- as.numeric(D[r, , drop = FALSE] %*% b)
  y <- spec$y[r]
  ll <- sum(spec$wt1[r] * w * (y * eta - log1p(exp(eta))))
  d <- b - mu_D[, j]
  ll - 0.5 * sum(d * (spec$P_D %*% d))
}

## Quadrature on a grid of +/- 8 Laplace SDs around the mode.
quad_moments <- function(j, w, mode, n_grid = 241L) {
  r <- spec$obs_rows[[j]]
  Xr <- D[r, , drop = FALSE]
  p <- stats::plogis(as.numeric(Xr %*% mode))
  H <- crossprod(Xr, Xr * (spec$wt1[r] * w * p * (1 - p))) + spec$P_D
  s <- sqrt(diag(solve(H)))
  g1 <- seq(mode[1] - 8 * s[1], mode[1] + 8 * s[1], length.out = n_grid)
  g2 <- seq(mode[2] - 8 * s[2], mode[2] + 8 * s[2], length.out = n_grid)
  lp <- outer(seq_along(g1), seq_along(g2), Vectorize(function(a, b) {
    log_target(c(g1[a], g2[b]), j, w)
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
  spec_w <- spec
  spec_w$wt1 <- spec$wt1 * (1 - w_rwt)
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
    q <- quad_moments(j, 1 - w_rwt, as.numeric(modes[a, ]))
    for (c in seq_len(ncol(D))) {
      x <- draws[, c, a]
      z <- (mean(x) - q$mean[c]) / sqrt(q$var[c] / n_draws)
      ratio <- stats::var(x) / q$var[c]
      z_var <- (ratio - 1) / sqrt(2 / (n_draws - 1))
      rows[[length(rows) + 1L]] <- data.frame(
        rwt = w_rwt, state = gl[j], which = names(pick)[a], coef = cn[c],
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

## Tempering check: for a data-rich state the posterior SD should grow
## roughly like 1 / sqrt(1 - rwt), so the SD ratio between the largest and
## smallest rwt is about sqrt((1 - rwt_lo) / (1 - rwt_hi)), a bit less
## because of the prior.
rich <- res[res$which == "rich", ]
rwt_lo <- min(rwt_values)
rwt_hi <- max(rwt_values)
for (cc in cn) {
  sd_lo <- rich$quad_sd[rich$rwt == rwt_lo & rich$coef == cc]
  sd_hi <- rich$quad_sd[rich$rwt == rwt_hi & rich$coef == cc]
  cat(sprintf("rich state, %-12s SD(rwt = %.2f) / SD(rwt = %.2f) = %.3f (likelihood only: %.3f)\n",
              cc, rwt_hi, rwt_lo, sd_hi / sd_lo,
              sqrt((1 - rwt_lo) / (1 - rwt_hi))))
}

z_all <- c(res$z_mean, res$z_var)
cat(sprintf(
  "\n%d z-scores: max |z| = %.2f, |z| > 3: %d (about %.2f expected by chance)\n",
  length(z_all), max(abs(z_all)), sum(abs(z_all) > 3),
  length(z_all) * 2 * stats::pnorm(-3)
))

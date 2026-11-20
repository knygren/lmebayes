## Prototype: predict the ridge TV distance TV(pi_rwt, pi_0) before sampling.
##
## pi_1 (the ridge at rwt = 1) replaces every group log-likelihood l_j by its
## second-order Taylor polynomial T_j at the mode. With fixed tau^2 and
## dNormal population priors it is an exact joint Gaussian over
## (b_1, ..., b_J, gamma), sampled here directly (no MCMC). With
##   h(b) = sum_j [l_j(b_j) - T_j(b_j)],
## pi_0 = pi_1 e^h / Z and pi_w = pi_1 e^((1-w) h) / Z_w, so the Laplace draws
## reweighted by e^h stand in for pi_0. Then, to first order in w,
##   TV(pi_w, pi_0) ~= w * kappa,   kappa = 0.5 * E_0 |h - E_0 h|,
## and rwt = tol / kappa keeps the ridge distance below tol.
##
## Printed: the check that the Gaussian's mean equals the ICM mode, weight
## diagnostics (ESS, Pareto k if 'loo' is installed), kappa with its batch
## spread, predicted TV at the rwt values already run, next to the TV the
## fitted draws gave, and the rwt each threshold would choose.
##
## Poisson, airbnb_small, same model as run_ridge_glmerb_poisson.R. Uses the
## installed lmebayesCore (no compilation). About 30 s.

suppressPackageStartupMessages(library(lmebayesCore))

n_lap <- 20000L
n_batch <- 5L
tol_values <- c(0.005, 0.01, 0.02)
rwt_check <- c(0.95, 0.75, 0.25, 0.05, 0.01)
## TV to rwt = 0 from importance reweighting of the 4000-chain fits
## (ridge_reweight_to_rwt0.R), same order as rwt_check.
tv_from_fits <- c(0.12, 0.0848, 0.0369, 0.00734, 0.00149)

data(airbnb_small, package = "bayesrules")
dat <- airbnb_small
dat$rating_c <- dat$rating - mean(dat$rating, na.rm = TRUE)
dat$walk_c <- dat$walk_score - mean(dat$walk_score, na.rm = TRUE)
dat <- dat[stats::complete.cases(dat[, c(
  "reviews", "rating_c", "walk_c", "neighborhood"
)]), ]
dat$neighborhood <- droplevels(factor(dat$neighborhood))
form <- reviews ~ walk_c + rating_c + (1 + rating_c || neighborhood)

design <- model_setup(form, data = dat, family = poisson())
ps <- Prior_Setup_GLMM(form, data = dat, family = poisson(), pop.pwt = 0.01)
pf <- pfamily_list(ps)
pt <- suppressWarnings(rglmerb(n = 1L, design = design, pfamily_list = pf,
                               family = poisson(), simulate = FALSE,
                               verbose = FALSE))

re_names <- design$groupef.names
group_levels <- levels(design$group)
J <- length(group_levels)
p_re <- length(re_names)
X_all <- as.matrix(design$D)[, re_names, drop = FALSE]
y_all <- as.numeric(design$y)
g_all <- as.character(design$group)
m_all <- as.matrix(pt$groupef.mode)[group_levels, re_names, drop = FALSE]

## Parameter layout: b_j (p_re each, group-major), then gamma_k per component
b_idx <- function(j) (j - 1L) * p_re + seq_len(p_re)
p_k <- vapply(re_names, function(k) ncol(as.matrix(design$W[[k]])), integer(1))
g_start <- J * p_re + c(0L, cumsum(p_k)[-length(p_k)])
g_idx <- lapply(seq_len(p_re), function(k) g_start[k] + seq_len(p_k[k]))
n_par <- J * p_re + sum(p_k)

Q <- matrix(0, n_par, n_par)
r <- numeric(n_par)
lik <- vector("list", J)

for (j in seq_len(J)) {
  rows <- g_all == group_levels[j]
  X <- X_all[rows, , drop = FALSE]
  y <- y_all[rows]
  m <- m_all[j, ]
  eta_m <- as.numeric(X %*% m)
  mu_m <- exp(eta_m)
  g <- as.numeric(crossprod(X, y - mu_m))
  C <- crossprod(X, X * mu_m)
  lik[[j]] <- list(X = X, y = y, m = m, g = g, C = C,
                   l_m = sum(y * eta_m - mu_m))
  ib <- b_idx(j)
  Q[ib, ib] <- Q[ib, ib] + C
  r[ib] <- r[ib] + as.numeric(C %*% m + g)
}

for (k in seq_len(p_re)) {
  nm <- re_names[k]
  Wk <- as.matrix(design$W[[nm]])[group_levels, , drop = FALSE]
  pl <- pf[[nm]]$prior_list
  tau2 <- as.numeric(pl$dispersion)
  for (j in seq_len(J)) {
    a <- numeric(n_par)
    a[b_idx(j)[k]] <- 1
    a[g_idx[[k]]] <- -Wk[j, ]
    Q <- Q + tcrossprod(a) / tau2
  }
  P_g <- solve(as.matrix(pl$Sigma))
  ig <- g_idx[[k]]
  Q[ig, ig] <- Q[ig, ig] + P_g
  r[ig] <- r[ig] + as.numeric(P_g %*% as.numeric(pl$mu))
}

R <- chol(Q)
theta_mean <- backsolve(R, forwardsolve(t(R), r))

mode_vec <- numeric(n_par)
for (j in seq_len(J)) mode_vec[b_idx(j)] <- m_all[j, ]
for (k in seq_len(p_re)) {
  mode_vec[g_idx[[k]]] <- as.numeric(pt$popef.mode[[re_names[k]]])
}
cat(sprintf("Laplace mean vs ICM mode: max |difference| = %.2e over %d parameters\n",
            max(abs(theta_mean - mode_vec)), n_par))

## Exact draws from pi_1 and h for each draw
Z <- matrix(stats::rnorm(n_par * n_lap), n_par, n_lap)
theta <- theta_mean + backsolve(R, Z)
h <- numeric(n_lap)
for (j in seq_len(J)) {
  L <- lik[[j]]
  B <- t(theta[b_idx(j), , drop = FALSE])
  E <- B %*% t(L$X)
  l_b <- as.numeric(E %*% L$y) - rowSums(exp(E)) - L$l_m
  D <- sweep(B, 2L, L$m)
  T_b <- as.numeric(D %*% L$g) - 0.5 * rowSums((D %*% L$C) * D)
  h <- h + (l_b - T_b)
}

norm_w <- function(lw) {
  w <- exp(lw - max(lw))
  w / sum(w)
}
kappa_hat <- function(hh) {
  W0 <- norm_w(hh)
  0.5 * sum(W0 * abs(hh - sum(W0 * hh)))
}

W0 <- norm_w(h)
ess0 <- 1 / sum(W0^2)
h_bar0 <- sum(W0 * h)
sd_h_pi1 <- stats::sd(h)
sd_h_pi0 <- sqrt(sum(W0 * (h - h_bar0)^2))
kappa <- kappa_hat(h)
batch <- split(seq_len(n_lap), rep(seq_len(n_batch), length.out = n_lap))
kappa_b <- vapply(batch, function(ix) kappa_hat(h[ix]), numeric(1))

cat(sprintf("\nLaplace draws: %d. ESS of e^h weights: %.0f (%.1f%%)\n",
            n_lap, ess0, 100 * ess0 / n_lap))
if (requireNamespace("loo", quietly = TRUE)) {
  k_hat <- suppressWarnings(loo::psis(h, r_eff = 1)$diagnostics$pareto_k)
  cat(sprintf("Pareto k of e^h weights: %.3f\n", k_hat))
}
cat(sprintf("sd(h): %.4f under pi_1, %.4f under pi_0 (reweighted)\n",
            sd_h_pi1, sd_h_pi0))
cat(sprintf("kappa = %.4f (batch SE %.4f); first-order rule 0.4 * sd_pi1(h) = %.4f\n",
            kappa, stats::sd(kappa_b) / sqrt(n_batch), 0.4 * sd_h_pi1))

tv_is <- vapply(rwt_check, function(w) {
  0.5 * sum(abs(W0 - norm_w((1 - w) * h)))
}, numeric(1))
cat("\n===== TV(pi_rwt, pi_0): predicted before sampling vs from the fits =====\n")
print(data.frame(
  rwt = rwt_check,
  kappa_x_rwt = signif(kappa * rwt_check, 3),
  laplace_IS = signif(tv_is, 3),
  rule_0.4sd = signif(0.4 * sd_h_pi1 * rwt_check, 3),
  from_fits = tv_from_fits
), row.names = FALSE)

cat("\n===== rwt chosen for each TV threshold (rwt = tol / kappa, capped at 0.95) =====\n")
print(data.frame(
  tol = tol_values,
  rwt = signif(pmin(tol_values / kappa, 0.95), 3)
), row.names = FALSE)

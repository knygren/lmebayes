## Where do the sampler's draws and the Laplace prediction disagree?
##
## For each saved fit at rwt = w, compares the spread of
## h_j(b) = l_j(b_j) - T_j(b_j) per group between the sampler's draws and
## the exact pi_w represented by Laplace draws reweighted by e^((1-w) h)
## (ridge_tv_predict_prototype.R). Also compares per-group SDs and
## 0.1% / 99.9% quantiles of the linear predictor at each group's mean
## rating, and kappa_w = 0.5 * E_w |h - E_w h| (the TV slope).
##
## Needs the fits saved by run_ridge_glmerb_poisson.R. About 1 min.

fits_path <- "C:/Rpackages/ridge_runs/ridge_fits_poisson.rds"
fits <- readRDS(fits_path)

proto <- new.env()
invisible(utils::capture.output(sys.source(
  "C:/Rpackages/lmebayes/data-raw/ridge_tv_predict_prototype.R", envir = proto
)))
re_names <- proto$re_names
group_levels <- proto$group_levels
J <- proto$J
listings <- vapply(proto$lik, function(L) nrow(L$X), integer(1))

h_group <- function(B, L) {
  E <- B %*% t(L$X)
  l_b <- as.numeric(E %*% L$y) - rowSums(exp(E)) - L$l_m
  D <- sweep(B, 2L, L$m)
  l_b - (as.numeric(D %*% L$g) - 0.5 * rowSums((D %*% L$C) * D))
}
group_draws <- function(fit, lv) {
  gdf <- fit$groupef
  gcol <- names(gdf)[2L]
  draws <- sort(unique(gdf$draw))
  sub <- gdf[as.character(gdf[[gcol]]) == lv, ]
  as.matrix(sub[match(draws, sub$draw), re_names])
}
wsd <- function(x, W) sqrt(sum(W * (x - sum(W * x))^2))
wq <- function(x, W, p) {
  o <- order(x)
  x[o][findInterval(p, cumsum(W[o])) + 1L]
}

hj_lap <- sapply(seq_len(J), function(j) {
  h_group(t(proto$theta[proto$b_idx(j), , drop = FALSE]), proto$lik[[j]])
})

for (nm in names(fits)) {
  w <- as.numeric(nm)
  fit <- fits[[nm]]
  m_fit <- as.matrix(fit$groupef.mode)[group_levels, re_names, drop = FALSE]
  Ww <- proto$norm_w((1 - w) * rowSums(hj_lap))

  hj_fit <- sapply(seq_len(J), function(j) {
    h_group(group_draws(fit, group_levels[j]), proto$lik[[j]])
  })
  h_fit <- rowSums(hj_fit)
  h_lap <- rowSums(hj_lap)
  kap_fit <- 0.5 * mean(abs(h_fit - mean(h_fit)))
  kap_lap <- 0.5 * sum(Ww * abs(h_lap - sum(Ww * h_lap)))

  cat(sprintf("\n========== rwt = %s (%d draws) ==========\n", nm, length(h_fit)))
  cat(sprintf("max |mode difference|, fit vs prototype: %.2e\n",
              max(abs(m_fit - proto$m_all))))
  cat(sprintf("kappa_w: fit %.4f, Laplace %.4f (ratio %.3f); sd(h): fit %.4f, Laplace %.4f\n",
              kap_fit, kap_lap, kap_fit / kap_lap, stats::sd(h_fit), wsd(h_lap, Ww)))

  tab <- do.call(rbind, lapply(seq_len(J), function(j) {
    L <- proto$lik[[j]]
    xbar <- colMeans(L$X)
    eta_fit <- as.numeric(group_draws(fit, group_levels[j]) %*% xbar)
    eta_lap <- as.numeric(t(proto$theta[proto$b_idx(j), , drop = FALSE]) %*% xbar)
    data.frame(
      group = group_levels[j], listings = listings[j],
      sd_hj_fit = stats::sd(hj_fit[, j]), sd_hj_lap = wsd(hj_lap[, j], Ww),
      sd_eta_fit = stats::sd(eta_fit), sd_eta_lap = wsd(eta_lap, Ww),
      q001_fit = stats::quantile(eta_fit, 0.001, names = FALSE),
      q001_lap = wq(eta_lap, Ww, 0.001),
      q999_fit = stats::quantile(eta_fit, 0.999, names = FALSE),
      q999_lap = wq(eta_lap, Ww, 0.999)
    )
  }))
  tab$ratio_hj <- tab$sd_hj_fit / tab$sd_hj_lap
  print(tab[order(-tab$listings), ], digits = 4, row.names = FALSE)
}

## How close are the ridge fits to rwt = 0? Importance reweighting.
##
## For the group effects the ridge target is an exponential tilt of the
## rwt = 0 posterior:
##   pi_w(b) = pi_0(b) * exp(-w * h(b)) / Z(w),   h(b) = sum_j [l_j(b_j) - T_j(b_j)],
## with T_j the second-order Taylor expansion of the group log-likelihood at
## the mode (the same g and C as .two_block_ridge_attach). The population
## block is unchanged, so draws from a fit at w reweighted by
## exp((w - w_target) * h(b)) target the rwt = w_target posterior.
## The weights also give the ESS and an estimate of TV(pi_w, pi_target).
##
## Run after run_ridge_glmerb_poisson.R in the same session: needs `fits`,
## `dat`, `skew`, `ex_kurt`. Poisson, no offset, unit weights. Assumes the
## rows of each popef matrix are in draw order.
## Validation: the 0.05 draws reweighted to 0.01 should match the 0.01 fit.

stopifnot(exists("fits"), exists("dat"))
re_cols <- c("(Intercept)", "rating_c")
rows_n <- sort(table(dat$neighborhood), decreasing = TRUE)
grp_levels <- names(rows_n)
y_all <- dat$reviews
X_all <- cbind(1, dat$rating_c)
g_all <- as.character(dat$neighborhood)

## Draws of one group, rows ordered by draw index
group_draws <- function(fit, lv) {
  gdf <- fit$groupef
  gcol <- names(gdf)[2L]
  draws <- sort(unique(gdf$draw))
  sub <- gdf[as.character(gdf[[gcol]]) == lv, ]
  as.matrix(sub[match(draws, sub$draw), re_cols])
}

## h(b) = sum_j [l_j(b_j) - T_j(b_j)] per draw (constants cancel)
ridge_gap <- function(fit) {
  m_all <- as.matrix(fit$groupef.mode)
  if (is.null(rownames(m_all))) stop("groupef.mode has no group rownames.")
  m_all <- m_all[, re_cols, drop = FALSE]
  out <- 0
  for (lv in rownames(m_all)) {
    r <- g_all == lv
    X <- X_all[r, , drop = FALSE]
    y <- y_all[r]
    B <- group_draws(fit, lv)
    m <- m_all[lv, ]
    eta_m <- as.numeric(X %*% m)
    mu_m <- exp(eta_m)
    g <- as.numeric(crossprod(X, y - mu_m))
    C <- crossprod(X, X * mu_m)
    E <- B %*% t(X)
    l_b <- as.numeric(E %*% y) - rowSums(exp(E)) - sum(y * eta_m - mu_m)
    D <- sweep(B, 2L, m)
    T_b <- as.numeric(D %*% g) - 0.5 * rowSums((D %*% C) * D)
    out <- out + (l_b - T_b)
  }
  out
}

## Normalized weights from rwt w_src to w_tgt, with diagnostics
reweight <- function(h, w_src, w_tgt) {
  lw <- (w_src - w_tgt) * h
  W <- exp(lw - max(lw))
  W <- W / sum(W)
  N <- length(W)
  list(W = W, ess = 1 / sum(W^2), tv = 0.5 * sum(abs(N * W - 1)) / N,
       sd_logw = stats::sd(lw))
}

wstats <- function(x, W) {
  mu <- sum(W * x)
  z <- x - mu
  v <- sum(W * z^2)
  c(mean = mu, sd = sqrt(v), skew = sum(W * z^3) / v^1.5,
    ex_kurt = sum(W * z^4) / v^2 - 3)
}

group_table <- function(fit, W) {
  do.call(rbind, lapply(grp_levels, function(lv) {
    B <- group_draws(fit, lv)
    do.call(rbind, lapply(re_cols, function(k) {
      data.frame(group = lv, listings = as.integer(rows_n[[lv]]), coef = k,
                 t(wstats(B[, k], W)), check.names = FALSE)
    }))
  }))
}

pop_table <- function(fit, W) {
  do.call(rbind, lapply(names(fit$popef), function(k) {
    M <- as.matrix(fit$popef[[k]])
    do.call(rbind, lapply(colnames(M), function(cn) {
      data.frame(component = k, coef = cn, t(wstats(M[, cn], W)))
    }))
  }))
}

w_src <- as.numeric(names(fits))
h_list <- lapply(fits, ridge_gap)

cat("===== distance to rwt = 0 (from importance weights) =====\n")
rw0 <- lapply(seq_along(fits), function(i) reweight(h_list[[i]], w_src[i], 0))
names(rw0) <- names(fits)
print(data.frame(
  rwt = w_src,
  N = vapply(h_list, length, integer(1)),
  ESS = round(vapply(rw0, `[[`, numeric(1), "ess")),
  TV_to_rwt0 = signif(vapply(rw0, `[[`, numeric(1), "tv"), 3),
  sd_log_weight = signif(vapply(rw0, `[[`, numeric(1), "sd_logw"), 3)
), row.names = FALSE)
cat("TV is roughly 0.4 * sd_log_weight when that is small; the estimate",
    "carries Monte Carlo noise, so treat tiny values as upper bounds.\n")

## rwt = 0 estimates from each fit, side by side
wide <- function(tabs, stat, keys) {
  out <- tabs[[1]][, keys]
  for (nm in names(tabs)) out[[paste0("from ", nm)]] <- tabs[[nm]][[stat]]
  out
}
grp_rw <- lapply(names(fits), function(nm) group_table(fits[[nm]], rw0[[nm]]$W))
names(grp_rw) <- names(fits)
pop_rw <- lapply(names(fits), function(nm) pop_table(fits[[nm]], rw0[[nm]]$W))
names(pop_rw) <- names(fits)

for (stat in c("mean", "sd", "skew")) {
  cat(sprintf("\n===== rwt = 0 estimate of population %s, from each fit =====\n", stat))
  print(wide(pop_rw, stat, c("component", "coef")), digits = 4, row.names = FALSE)
}
for (stat in c("mean", "sd", "skew")) {
  cat(sprintf("\n===== rwt = 0 estimate of group %s, from each fit =====\n", stat))
  print(wide(grp_rw, stat, c("group", "listings", "coef")), digits = 4,
        row.names = FALSE)
}

## Linear extrapolation in rwt from the raw (unweighted) 0.05 and 0.01 means
if (all(c("0.05", "0.01") %in% names(fits))) {
  N05 <- length(h_list[["0.05"]])
  N01 <- length(h_list[["0.01"]])
  raw05 <- group_table(fits[["0.05"]], rep(1 / N05, N05))
  raw01 <- group_table(fits[["0.01"]], rep(1 / N01, N01))
  cat("\n===== group means: raw 0.01, linear extrapolation (5*E01 - E05)/4,",
      "and reweighted from 0.01 =====\n")
  print(data.frame(
    raw05[, c("group", "listings", "coef")],
    raw_0.01 = raw01$mean,
    extrap_0 = (5 * raw01$mean - raw05$mean) / 4,
    reweight_0 = grp_rw[["0.01"]]$mean
  ), digits = 5, row.names = FALSE)

  ## Validation: 0.05 draws reweighted to 0.01 against the 0.01 fit
  rw_val <- reweight(h_list[["0.05"]], 0.05, 0.01)
  val <- group_table(fits[["0.05"]], rw_val$W)
  se01 <- raw01$sd / sqrt(N01)
  cat(sprintf(paste0("\n===== validation: 0.05 draws reweighted to 0.01 ",
                     "(ESS %.0f of %d) vs the 0.01 fit =====\n"),
              rw_val$ess, N05))
  print(data.frame(
    raw01[, c("group", "listings", "coef")],
    fit_0.01 = raw01$mean,
    reweighted_0.05 = val$mean,
    z = (val$mean - raw01$mean) / (se01 * sqrt(1 + N01 / rw_val$ess))
  ), digits = 4, row.names = FALSE)
}

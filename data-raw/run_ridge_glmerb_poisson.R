## Full two-block Gibbs runs of the ridge model, Poisson.
##
## rglmerb() on airbnb_small, reviews ~ walk_c + rating_c +
## (1 + rating_c || neighborhood), n independent chains, one fit per rwt.
## Each fit runs the pilot stage and then the main stage for the
## Corollary 5.1 sweep count at that rwt (more sweeps for smaller rwt).
## The population priors are dNormal, so the compiled sweep
## (draw_engine "rGLMM_sweep_cpp") runs, with the groups likelihood raised
## to 1 - rwt and the ridge term in each group's prior.
##
## Reported per rwt: wall time, sweep count, ridge rate, population
## coefficient means, SDs, skewness and excess kurtosis, and the same for the
## group coefficients of every neighborhood in the model (sorted by
## listings, one column per rwt). Fits are kept in `fits`, named by rwt, and
## saved to `fits_path`. The ridge model's bias
## shrinks with rwt, so the summaries should settle as rwt decreases.
## glmer estimates are printed as a reference point.
##
## Manual only (long for small rwt). Run from anywhere:
##   Rscript C:/Rpackages/lmebayes/data-raw/run_ridge_glmerb_poisson.R
##   Rscript C:/Rpackages/lmebayes/data-raw/run_ridge_glmerb_poisson.R 0.1 0.05
## Arguments are rwt values (default 0.95, 0.75, 0.25, 0.05, 0.01). n is set below.

suppressPackageStartupMessages(library(pkgload))
pkgload::load_all("C:/Rpackages/lmebayesCore", quiet = TRUE)

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This script requires the 'bayesrules' package.", call. = FALSE)
}

n_chains <- 4000L
## Outside both package repositories; read by ridge_h_spread_compare.R
fits_path <- "C:/Rpackages/ridge_runs/ridge_fits_poisson.rds"
args <- commandArgs(trailingOnly = TRUE)
rwt_values <- if (length(args)) as.numeric(args) else c(0.95, 0.75, 0.25, 0.05, 0.01)
if (anyNA(rwt_values) || any(rwt_values <= 0 | rwt_values >= 1)) {
  stop("rwt values must be in (0, 1).", call. = FALSE)
}

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

rows_n <- sort(table(dat$neighborhood), decreasing = TRUE)
grp_levels <- names(rows_n)
cat(sprintf("Model: %d listings in %d neighborhoods (all are sampled and reported).\n",
            nrow(dat), length(grp_levels)))
print(data.frame(neighborhood = grp_levels, listings = as.integer(rows_n)),
      row.names = FALSE)

skew <- function(x) {
  z <- x - mean(x)
  mean(z^3) / mean(z^2)^1.5
}
ex_kurt <- function(x) {
  z <- x - mean(x)
  mean(z^4) / mean(z^2)^2 - 3
}

popef_summary <- function(fit) {
  do.call(rbind, lapply(names(fit$popef), function(k) {
    M <- as.matrix(fit$popef[[k]])
    data.frame(component = k, coef = colnames(M),
               mean = colMeans(M), sd = apply(M, 2, stats::sd),
               skew = apply(M, 2, skew), ex_kurt = apply(M, 2, ex_kurt),
               mode = as.numeric(fit$popef.mode[[k]]),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

groupef_summary <- function(fit) {
  g <- fit$groupef
  gcol <- names(g)[2L]
  coefs <- setdiff(names(g), c("draw", gcol))
  do.call(rbind, lapply(grp_levels, function(lv) {
    sub <- g[as.character(g[[gcol]]) == lv, coefs, drop = FALSE]
    data.frame(group = lv, listings = as.integer(rows_n[[lv]]), coef = coefs,
               mean = colMeans(sub), sd = apply(sub, 2, stats::sd),
               skew = apply(sub, 2, skew), ex_kurt = apply(sub, 2, ex_kurt),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

fits <- list()
runs <- list()
pop_rows <- list()
grp_rows <- list()
for (w in rwt_values) {
  cat(sprintf("\n===== rwt = %g, n = %d chains =====\n", w, n_chains))
  t0 <- proc.time()[["elapsed"]]
  fit <- lmebayesCore::rglmerb(
    n            = n_chains,
    design       = design,
    pfamily_list = pf,
    family       = poisson(),
    rwt          = w,
    verbose      = FALSE,
    progbar      = TRUE
  )
  secs <- proc.time()[["elapsed"]] - t0
  fits[[format(w)]] <- fit
  ci <- fit$convergence_info
  engine <- ci$draw_engine
  if (!identical(engine, "rGLMM_sweep_cpp")) {
    warning("draw_engine is ", engine, "; the ridge prior is not applied there.",
            call. = FALSE)
  }
  runs[[length(runs) + 1L]] <- data.frame(
    rwt = w,
    seconds = round(secs, 1),
    pilot_n = if (is.null(fit$pilot$n)) 0L else as.integer(fit$pilot$n),
    m_convergence = as.integer(fit$m_convergence),
    ridge_m = if (is.null(ci$ridge)) NA_integer_ else as.integer(ci$ridge$m_convergence),
    ridge_rate = if (is.null(ci$ridge)) NA_real_ else ci$ridge$lambda_star,
    engine = engine,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %.1f s, main sweeps = %d (ridge rate %.5f)\n", secs,
              as.integer(fit$m_convergence),
              if (is.null(ci$ridge)) NA_real_ else ci$ridge$lambda_star))
  pop_rows[[length(pop_rows) + 1L]] <- cbind(rwt = w, popef_summary(fit))
  grp_rows[[length(grp_rows) + 1L]] <- cbind(rwt = w, groupef_summary(fit))
}

dir.create(dirname(fits_path), showWarnings = FALSE, recursive = TRUE)
saveRDS(fits, fits_path)
cat(sprintf("\nFits saved to %s\n", fits_path))

cat("\n===== runs =====\n")
print(do.call(rbind, runs), row.names = FALSE)

cat("\n===== population coefficients (mean, sd over chains; mode = ICM) =====\n")
pop <- do.call(rbind, pop_rows)
print(pop[order(pop$component, pop$coef, -pop$rwt), ], digits = 4,
      row.names = FALSE)

cat(sprintf(paste0(
  "\nMonte Carlo SE under normality with n = %d: skew %.3f, ex_kurt %.3f\n"),
  n_chains, sqrt(6 / n_chains), sqrt(24 / n_chains)))

## Group tables: one row per group and coefficient (groups sorted by
## listings), one column per rwt.
grp <- do.call(rbind, grp_rows)
grp_wide <- function(stat) {
  keys <- unique(grp[, c("group", "listings", "coef")])
  keys <- keys[order(match(keys$group, grp_levels), keys$coef), ]
  for (w in rwt_values) {
    s <- grp[grp$rwt == w, ]
    keys[[paste0("rwt=", format(w))]] <-
      s[[stat]][match(paste(keys$group, keys$coef), paste(s$group, s$coef))]
  }
  keys
}
gm <- if (requireNamespace("lme4", quietly = TRUE)) {
  lme4::glmer(form, data = dat, family = poisson())
}
for (stat in c("mean", "sd", "skew", "ex_kurt")) {
  cat(sprintf("\n===== group coefficients: %s =====\n", stat))
  tab <- grp_wide(stat)
  if (identical(stat, "mean") && !is.null(gm)) {
    co_gm <- stats::coef(gm)$neighborhood
    tab$glmer <- mapply(function(g, k) co_gm[g, k], tab$group, tab$coef)
  }
  print(tab, digits = 4, row.names = FALSE)
}

if (!is.null(gm)) {
  cat("\n===== glmer fixed effects =====\n")
  print(lme4::fixef(gm), digits = 4)
}

## Demo: lmerb() full workflow on bayesrules::big_word_club
##
## Replica of the original ?lmerb example before it was simplified (the
## man-page example now fits a smaller intercept + free_reduced_lunch +
## distracted_ppvt model).  Preserved as a demo because the full run (1000
## draws plus the factor-level diagnostics) takes on the order of a minute.
## This demo keeps the complete workflow: model_setup(),
## Prior_Setup_lmebayes(), pfamily_list(), lmerb(), then draws-vs-ICM
## z-tests and lmer/mu_all/lmerb factor-level comparisons.
##
##   demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")

for (pkg in c("bayesrules", "coda")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("This demo requires the '%s' package.", pkg), call. = FALSE)
  }
}

## lmerb: school random effects on bayesrules::big_word_club
##
## Compare factor-level means: lmer coef(), Block 1 prior mean mu_all,
## and posterior means of lmerb draws (averaged over draw within each level).
##
## Student-level random slopes: distracted_a1 and distracted_ppvt (distraction
## during assessment 1 vs during PPVT). Both have fixed main effects + school
## random slopes. Cross-level moderation: free_reduced_lunch:distracted_a1
## (school lunch status moderates the distracted_a1 random slope).
##
## Workflow: model_setup(), Prior_Setup_lmebayes(), pfamily_list(), then
## lmerb(pfamily_list = , dispersion_ranef = ).

data(big_word_club, package = "bayesrules")

dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c(
      "score_ppvt", "distracted_a1", "distracted_ppvt",
      "private_school", "title1", "free_reduced_lunch", "school_id"
    )])
)


form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_ppvt + distracted_a1 +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)



design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
cat("\n=== Prior_Setup_lmebayes ===\n\n")
print(ps)

## Defaults: tv_tol = 0.01 (each stored draw within 0.01 TV of the exact
## joint posterior; m_convergence derived from the Nygren (2020) Theorem 3
## bound), no fixed seed.
fit <- lmerb(
  form_lmer,
  data = dat,
  pfamily_list = pfamily_list(ps),
  dispersion_ranef = ps$dispersion_ranef,
  n = 1000L
)

cat("\n=== summary(fit) ===\n\n")
print(summary(fit))

grp_col  <- fit$model_setup$group_name
re_names <- fit$model_setup$re_coef_names
grp_levs <- rownames(coef(fit$lmer)[[grp_col]])

## --- Posterior means of fixed effects from Block 2 draws -----------------
## fit$fixef.means[[k]]: MCMC mean of the n Block 2 gamma_k draws.
## Directly comparable to lme4::fixef() (same parameter scale).
## fit$fixef.mode: exact ICM posterior mean (used as H0 for z-test below).
## z = (draws mean - ICM mean) / (draws SD / sqrt(n)); flag |z| > 2 with *.
cat("\n=== Posterior means of fixed effects (from Block 2 draws) ===\n\n")
n_draws     <- nrow(fit$fixef[[re_names[1L]]])
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            "RE component", "parameter", "draws mean", "draws SD", "SE(mean)", "ICM mean", "z"))
cat(sprintf("  %-18s  %-28s  %10s  %10s  %10s  %10s  %7s\n",
            strrep("-", 18L), strrep("-", 28L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L),
            strrep("-", 10L), strrep("-", 7L)))
for (k in re_names) {
  dm_k  <- fit$fixef.means[[k]]
  sd_k  <- apply(fit$fixef[[k]], 2L, sd)
  se_k  <- sd_k / sqrt(n_draws)
  icm_k <- fit$fixef.mode[[k]]
  for (nm in names(dm_k)) {
    z_val <- (dm_k[[nm]] - icm_k[[nm]]) / se_k[[nm]]
    flag  <- if (abs(z_val) > 2) " *" else "  "
    cat(sprintf("  %-18s  %-28s  %10.4f  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                k, nm, dm_k[[nm]], sd_k[[nm]], se_k[[nm]], icm_k[[nm]], z_val, flag))
  }
}
cat("  (* |z| > 2: draws mean inconsistent with exact ICM posterior mean)\n")
cat("\n")

## --- Z-test: MCMC mean vs ICM posterior mean for random effects -------------
## fit$ranef.mode: exact ICM posterior mean (J x p_re), rows = group levels.
## MCMC mean computed as the per-group average of the n Block 1 draws.
## z = (MCMC mean - ICM mean) / (draws SD / sqrt(n))
## With J * p_re tests we expect a few |z| > 2 by chance; flag |z| > 3.

cat("=== Random effects: MCMC mean vs ICM posterior mean ===\n\n")

## Compute per-group MCMC means and SDs from the stored coefficients.
## simplify=FALSE forces a named list regardless of output length.
re_draws_mean <- tapply(
  seq_len(nrow(fit$coefficients)),
  fit$coefficients[[grp_col]],
  function(idx) colMeans(fit$coefficients[idx, re_names, drop = FALSE]),
  simplify = FALSE
)
re_draws_sd <- tapply(
  seq_len(nrow(fit$coefficients)),
  fit$coefficients[[grp_col]],
  function(idx) apply(fit$coefficients[idx, re_names, drop = FALSE], 2L, sd),
  simplify = FALSE
)
n_re_draws <- n_draws   # same n as fixed effects

cat(sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
            "group", "RE component", "MCMC mean", "ICM mean", "SE(mean)", "z"))
cat(sprintf("  %-6s  %-18s  %10s  %10s  %10s  %6s\n",
            strrep("-", 6L), strrep("-", 18L),
            strrep("-", 10L), strrep("-", 10L), strrep("-", 10L), strrep("-", 6L)))

icm_b <- fit$ranef.mode   # J x p_re, rownames = group levels

n_flagged <- 0L
for (lev in grp_levs) {
  lev_chr <- as.character(lev)
  for (k in re_names) {
    mcmc_m <- re_draws_mean[[lev_chr]][[k]]
    mcmc_s <- re_draws_sd[[lev_chr]][[k]]
    icm_m  <- icm_b[lev_chr, k]
    se_val <- mcmc_s / sqrt(n_re_draws)
    z_val  <- (mcmc_m - icm_m) / se_val
    flag   <- if (abs(z_val) > 3) " *" else "  "
    if (abs(z_val) > 3) n_flagged <- n_flagged + 1L
    cat(sprintf("  %-6s  %-18s  %10.4f  %10.4f  %10.4f  %6.2f%s\n",
                lev_chr, k, mcmc_m, icm_m, se_val, z_val, flag))
  }
}
total_tests <- length(grp_levs) * length(re_names)
cat(sprintf(
  "\n  %d of %d tests flagged |z| > 3  (expected ~%.1f by chance at 0.3%% level)\n",
  n_flagged, total_tests, total_tests * 0.003
))
cat("  (* |z| > 3: MCMC mean inconsistent with exact ICM posterior mean)\n\n")

rename_re_cols <- function(df, re_names, suffix) {
  idx <- match(re_names, names(df))
  if (anyNA(idx)) {
    stop("missing RE columns: ", paste(re_names[is.na(idx)], collapse = ", "))
  }
  names(df)[idx] <- paste0(re_names, suffix)
  df
}

fe_lmer <- lme4::fixef(fit$lmer)
coef_raw_df <- as.data.frame(coef(fit$lmer)[[grp_col]])
## anchor[k]: fixef(lmer) for terms with a population fixed effect.
coef_anchor <- vapply(re_names, function(k) {
  if (k == "(Intercept)") {
    unname(fe_lmer["(Intercept)"])
  } else if (k %in% names(fe_lmer)) {
    unname(fe_lmer[k])
  } else {
    unname(fit$fixef.mode[[k]]["(Intercept)"])
  }
}, numeric(1L))

## lmer: total random-effects coefficients at each grouping-factor level
lmer_by_level <- coef_raw_df
lmer_by_level[[grp_col]] <- factor(rownames(lmer_by_level), levels = grp_levs)
rownames(lmer_by_level) <- NULL
lmer_by_level <- rename_re_cols(lmer_by_level, re_names, "_lmer")

mu_mat <- as.matrix(fit$fixef.mu)
mu_by_level <- as.data.frame(t(mu_mat), stringsAsFactors = FALSE)
mu_by_level[[grp_col]] <- factor(rownames(mu_by_level), levels = grp_levs)
rownames(mu_by_level) <- NULL
mu_by_level <- rename_re_cols(mu_by_level, re_names, "_mu_all")

## lmer aligned to lmebayes: mu_all + (coef - anchor)
lmer_full_df <- coef_raw_df
for (k in re_names) {
  mu_k <- mu_mat[k, rownames(coef_raw_df), drop = TRUE]
  lmer_full_df[[k]] <- mu_k + (coef_raw_df[[k]] - coef_anchor[k])
}
lmer_full_by_level <- lmer_full_df
lmer_full_by_level[[grp_col]] <- factor(rownames(lmer_full_df), levels = grp_levs)
rownames(lmer_full_by_level) <- NULL
lmer_full_by_level <- rename_re_cols(lmer_full_by_level, re_names, "_lmer_full")

## lmerb: posterior mean at each grouping-factor level (mean over draw)
lmerb_by_level <- aggregate(
  fit$coefficients[, re_names, drop = FALSE],
  by = list(fit$coefficients[[grp_col]]),
  FUN = mean,
  simplify = TRUE
)
names(lmerb_by_level)[1L] <- grp_col
lmerb_by_level[[grp_col]] <- factor(lmerb_by_level[[grp_col]], levels = grp_levs)
lmerb_by_level <- rename_re_cols(lmerb_by_level, re_names, "_lmerb")

## Factor-level comparison (one row per school)
level_means <- merge(lmer_by_level, lmer_full_by_level, by = grp_col, sort = FALSE)
level_means <- merge(level_means, mu_by_level, by = grp_col, sort = FALSE)
level_means <- merge(level_means, lmerb_by_level, by = grp_col, sort = FALSE)
level_means <- level_means[order(level_means[[grp_col]]), , drop = FALSE]

cat("\nMean coefficient across schools (raw coef vs mu_all vs lmer_full vs lmerb):\n")
avg_row <- data.frame(
  term = re_names,
  lmer_raw = vapply(re_names, function(nm) {
    mean(level_means[[paste0(nm, "_lmer")]])
  }, numeric(1L)),
  mu_all = vapply(re_names, function(nm) {
    mean(level_means[[paste0(nm, "_mu_all")]])
  }, numeric(1L)),
  lmer_full = vapply(re_names, function(nm) {
    mean(level_means[[paste0(nm, "_lmer_full")]])
  }, numeric(1L)),
  lmerb = vapply(re_names, function(nm) {
    mean(level_means[[paste0(nm, "_lmerb")]])
  }, numeric(1L)),
  row.names = NULL
)
avg_row[-1L] <- lapply(avg_row[-1L], function(x) round(x, 3))
print(avg_row)

cat("\nFactor-level means: lmer (raw coef), lmer_full, mu_all, lmerb:\n")
show_cols <- c(
  grp_col,
  paste0(re_names, "_lmer"),
  paste0(re_names, "_lmer_full"),
  paste0(re_names, "_mu_all"),
  paste0(re_names, "_lmerb")
)
out_means <- level_means[, show_cols, drop = FALSE]
num_cols <- setdiff(show_cols, grp_col)
out_means[num_cols] <- lapply(out_means[num_cols], function(x) round(as.numeric(x), 3))
print(out_means)

## Long format: one row per (factor level, RE term)
## Compare lmerb to lmer_full (lmebayes-aligned lmer), not raw coef(lmer).
level_long <- do.call(rbind, lapply(re_names, function(nm) {
  lmer_full_v <- level_means[[paste0(nm, "_lmer_full")]]
  mu_all_v    <- level_means[[paste0(nm, "_mu_all")]]
  lmerb_v     <- level_means[[paste0(nm, "_lmerb")]]
  data.frame(
    level      = level_means[[grp_col]],
    term       = nm,
    lmer_raw   = level_means[[paste0(nm, "_lmer")]],
    lmer_full  = lmer_full_v,
    mu_all     = mu_all_v,
    lmerb      = lmerb_v,
    u_lmer     = lmer_full_v - mu_all_v,
    u_lmerb    = lmerb_v - mu_all_v,
    diff_bf    = lmerb_v - lmer_full_v,
    diff_u     = (lmerb_v - mu_all_v) - (lmer_full_v - mu_all_v),
    stringsAsFactors = FALSE
  )
}))

cat("\nFactor-level comparison by term (all schools):\n")
cat("  lmer_raw   = coef(lmer)\n")
cat("  lmer_full  = mu_all + (lmer_raw - anchor)  [full b_j in lmebayes notation]\n")
cat("  u_lmer     = lmer_full - mu_all  (school RE deviation from prior mean)\n")
cat("  u_lmerb    = lmerb - mu_all      (same deviation from posterior mean)\n")
cat("  diff_bf    = lmerb - lmer_full;  diff_u = u_lmerb - u_lmer (= diff_bf)\n")
cat("  distracted_ppvt / distracted_a1: fixed + RE; anchor = fixef()\n")
cat("  distracted_a1 also moderated by free_reduced_lunch (mu_all varies by school)\n\n")
out_long <- level_long
num_long <- c(
  "lmer_raw", "lmer_full", "mu_all", "lmerb",
  "u_lmer", "u_lmerb", "diff_bf", "diff_u"
)
out_long[num_long] <- lapply(out_long[num_long], function(x) round(x, 3))
print(out_long)

cat("\nMean |lmerb - lmer_full| across factor levels, by term:\n")
print(round(tapply(abs(level_long$diff_bf), level_long$term, mean), 3))
cat("\nMean |u_lmerb - u_lmer| (= |diff_bf|) across factor levels, by term:\n")
print(round(tapply(abs(level_long$diff_u), level_long$term, mean), 3))
cat("\nMean |lmerb - lmer_raw| across factor levels, by term:\n")
print(round(tapply(
  abs(level_long$lmerb - level_long$lmer_raw),
  level_long$term,
  mean
), 3))

cat("\n--- distracted_ppvt: fixed + RE (replaces age_c) ---\n")
cat(sprintf(
  "  anchor = fixef distracted_ppvt: %.3f; mu_all (constant): %.3f\n",
  coef_anchor["distracted_ppvt"], mu_mat["distracted_ppvt", 1L]
))
ppvt_long <- subset(level_long, term == "distracted_ppvt")
cat("  Cor(lmer_full, lmerb): ", round(cor(ppvt_long$lmer_full, ppvt_long$lmerb), 3), "\n")
cat("  Cor(u_lmer, u_lmerb):   ", round(cor(ppvt_long$u_lmer, ppvt_long$u_lmerb), 3), "\n")
cat("  Mean |diff_bf|:         ", round(mean(abs(ppvt_long$diff_bf)), 3), "\n")

## coda: one chain per factor level (10000 draws x p_re per school)
coef_split <- split(
  fit$coefficients,
  factor(fit$coefficients[[grp_col]], levels = grp_levs)
)
mcmc_by_level <- coda::mcmc.list(lapply(coef_split, function(d) {
  coda::mcmc(as.matrix(d[, re_names, drop = FALSE]))
}))
ex_level <- as.character(grp_levs[1L])
cat("\ncoda summary for factor level", ex_level, ":\n")
print(summary(mcmc_by_level[[ex_level]]))

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
## Workflow: model_setup(), then Prior_Setup_lmebayes(), then lmerb().

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
  distracted_a1 + distracted_ppvt +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup ===\n\n")
print(design)

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)



fit <- lmerb(
  form_lmer,
  data = dat,
  measurement_prior_list = ps,
  n = 10L,
  seed = 42L
)

grp_col  <- fit$model_setup$group_name
re_names <- fit$model_setup$re_coef_names
grp_levs <- rownames(coef(fit$lmer)[[grp_col]])

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
    unname(fit$fixef[[k]]["(Intercept)"])
  }
}, numeric(1L))

## lmer: total random-effects coefficients at each grouping-factor level
lmer_by_level <- coef_raw_df
lmer_by_level[[grp_col]] <- factor(rownames(lmer_by_level), levels = grp_levs)
rownames(lmer_by_level) <- NULL
lmer_by_level <- rename_re_cols(lmer_by_level, re_names, "_lmer")

mu_mat <- as.matrix(fit$mu_all)
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

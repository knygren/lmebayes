## Demo: lmerb() with dGamma() sigma^2 and ING Block~2 on big_word_club
##
## Gaussian LMM route rLMMindepNormalGamma_reg_estimated_vcov(): random sigma^2
## (dGamma dispersion_ranef) plus sampled tau^2_k (dIndependent_Normal_Gamma
## Block~2). Combines demo/Ex_21_lmerb_ING_BigWordClub.R with Ex_23 case 4.
##
## TEMP: subsets to algebraically full-rank schools (Block~1 ING per-group
## requires full-rank Z_j). Also excludes school_id 2 and 18 (envelope sign
## violation in per-group ING). Remove when shared sigma^2 Block~1 is in place.
##
##   demo("Ex_25_lmerb_dGamma_ING_BigWordClub", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This demo requires the 'lme4' package.", call. = FALSE)
}

## TEMP: trace Block~1 per-group ING progress (remove when debugging done)
options(glmbayesCore.debug_block1_ing_levels = TRUE)

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

design_all <- model_setup(form_lmer, data = dat)
full_rank_schools <- names(design_all$re_rank)[design_all$re_rank]
cat(sprintf(
  "\n=== Full-rank filter: %d of %d schools kept ===\n",
  length(full_rank_schools),
  length(design_all$re_rank)
))
if (length(full_rank_schools) < length(design_all$re_rank)) {
  cat(
    "  Dropped:",
    paste(names(design_all$re_rank)[!design_all$re_rank], collapse = ", "),
    "\n"
  )
}
dat <- subset(dat, school_id %in% full_rank_schools)
dat$school_id <- droplevels(dat$school_id)

## TEMP: schools 2 and 18 trigger ING envelope sign violation (UB2 < 0); drop and retest
temp_drop_schools <- c("18", "2")
drop <- intersect(temp_drop_schools, levels(dat$school_id))
if (length(drop)) {
  cat(sprintf(
    "\n=== TEMP: excluding school_id %s (Block~1 ING envelope failure) ===\n",
    paste(drop, collapse = ", ")
  ))
  dat <- subset(dat, !as.character(school_id) %in% drop)
  dat$school_id <- droplevels(dat$school_id)
}

design <- model_setup(form_lmer, data = dat)
cat("\n=== model_setup (full-rank schools only) ===\n\n")
print(design)
stopifnot(all(design$re_rank))

ps <- Prior_Setup_lmebayes(
  form_lmer,
  data           = dat,
  pwt            = 0.01,
  pwt_dispersion = 0.2
)
cat("\n=== Prior_Setup_lmebayes (ING calibration) ===\n\n")
print(ps)

pf <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")

m_disp <- ps$ing_prior_measurement
disp_pf <- dGamma(
  shape          = m_disp$shape,
  rate           = m_disp$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m_disp$disp_lower,
  disp_upper     = m_disp$disp_upper
)

cat("\n=== lmer reference fit ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
print(summary(fit_lmer))

fit <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = disp_pf,
  n                = 3000L,
  gap_tol          = 0.05,
  mode_gap_max     = 1.0,
  diag_sweeps      = FALSE
)

stopifnot(identical(fit$prior$dispersion_mode, "gamma"))
stopifnot(isTRUE(fit$prior$any_non_normal))
stopifnot(!is.null(fit$fixef.dispersion))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(fit$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))
stopifnot(is.finite(fit$pilot_chisq$p_value))
stopifnot(!is.null(fit$sweep_history$pilot))
stopifnot(!is.null(fit$sweep_history$main))

cat("\n=== summary(lmerb fit) ===\n\n")
print(summary(fit))

cat("\n=== Block~2 sweep summaries (pilot, then main) ===\n\n")
print(fit$sweep_history$pilot)
print(fit$sweep_history$main)

re_names <- fit$model_setup$re_coef_names

for (k in re_names) {
  pr_k <- pf[[k]]$prior_list
  t2   <- fit$fixef.dispersion[, k]
  stopifnot(
    all(is.finite(t2)), all(t2 > 0),
    all(t2 >= pr_k$disp_lower),
    all(t2 <= pr_k$disp_upper),
    stats::sd(t2) > 0
  )
  cat(sprintf(
    "\n%s tau^2: post mean = %.4f  [window (%.4f, %.4f)]\n",
    k,
    fit$fixef.dispersion.mean[[k]],
    pr_k$disp_lower,
    pr_k$disp_upper
  ))
}

cat(sprintf(
  "\nPilot vs plug-in start (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  fit$convergence$m_convergence_pilot
))

cn <- unlist(lapply(re_names, function(k) {
  paste0(k, "::", colnames(fit$fixef[[k]]))
}))
beta_bar <- unlist(lapply(re_names, function(k) fit$fixef.means[[k]]))
theta_plug <- unlist(lapply(re_names, function(k) fit$fixef.mode[[k]]))
theta_pilot <- unlist(lapply(re_names, function(k) {
  nms <- colnames(fit$fixef[[k]])
  unname(fit$fixef.init[[k]][nms])
}))
names(beta_bar) <- names(theta_plug) <- names(theta_pilot) <- cn

block2_cmp <- data.frame(
  plug_in     = unname(theta_plug),
  pilot_mean  = unname(theta_pilot),
  mcmc_mean   = unname(beta_bar),
  row.names   = cn,
  check.names = FALSE
)
cat("\n=== Block 2 hyperparameters (plug-in / pilot / MCMC) ===\n\n")
print(round(block2_cmp, 4))

cat("\n=== Random effects: lmer reference vs lmerb chain mean ===\n\n")
lmebayes:::print_mer_bayes_re_compare(fit)

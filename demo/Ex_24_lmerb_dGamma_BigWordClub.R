## Demo: lmerb() with dGamma() observation dispersion on big_word_club
##
## Gaussian LMM route rLMMindepNormalGamma_reg_known_vcov(): random sigma^2
## (dGamma dispersion_ranef) with fixed Block~2 tau^2_k (all dNormal pfamily).
## Same model as demo/Ex_12_lmerb_BigWordClub.R; compare Ex_23 case 3
## (simulate = FALSE joint mode only).
##
## TEMP: subsets to algebraically full-rank schools (Block~1 ING per-group
## requires full-rank Z_j). Also excludes school_id 2 and 18 (envelope sign
## violation in per-group ING). Remove when shared sigma^2 Block~1 is in place.
##
##   demo("Ex_24_lmerb_dGamma_BigWordClub", package = "lmebayes")

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

## TEMP: school 18 triggers ING envelope sign violation (UB2 < 0); drop and retest
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

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
cat("\n=== Prior_Setup_lmebayes ===\n\n")
print(ps)

pf <- pfamily_list(ps)

m_disp <- ps$ing_prior_measurement
disp_pf <- dGamma(
  shape          = m_disp$shape,
  rate           = m_disp$rate,
  beta           = matrix(0, 1, 1, dimnames = list("(Intercept)", NULL)),
  Inv_Dispersion = TRUE,
  disp_lower     = m_disp$disp_lower+0.25*(m_disp$disp_upper-m_disp$disp_lower),
  disp_upper     = m_disp$disp_upper-0.25*(m_disp$disp_upper-m_disp$disp_lower)
)

cat("\n=== lmer reference fit ===\n\n")
fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
print(summary(fit_lmer))

cat(sprintf(
  "\n=== dGamma sigma^2 prior mean (rate/(shape-1)): %.4f (REML sigma^2: %.4f) ===\n\n",
  m_disp$rate / (m_disp$shape - 1),
  stats::sigma(fit_lmer)^2
))

fit <- lmerb(
  form_lmer,
  data             = dat,
  pfamily_list     = pf,
  dispersion_ranef = disp_pf,
  n                = 1000L
)

stopifnot(identical(fit$prior$dispersion_mode, "gamma"))
stopifnot(!isTRUE(fit$prior$any_non_normal))
stopifnot(is.null(fit$fixef.dispersion))
stopifnot(!is.null(fit$pilot_chisq))
stopifnot(fit$pilot_chisq$n_pilot > 0L)
stopifnot(identical(fit$pilot_chisq$n_pilot, fit$convergence$n_pilot))
stopifnot(is.finite(fit$pilot_chisq$p_value))
stopifnot(!is.null(fit$sweep_history$main))

cat("\n=== summary(lmerb fit) ===\n\n")
print(summary(fit))

cat(sprintf(
  "\nPilot vs mode (chi-squared): p = %.4g (n_pilot = %d, m_convergence_pilot = %d)\n",
  fit$pilot_chisq$p_value,
  fit$pilot_chisq$n_pilot,
  fit$convergence$m_convergence_pilot
))

re_names <- fit$model_setup$re_coef_names
n_draws  <- nrow(fit$fixef[[re_names[1L]]])
stopifnot(identical(n_draws, 1000L))

cat("\n=== Block 2 fixed effects: draws mean vs ICM mean ===\n\n")
for (k in re_names) {
  dm_k  <- fit$fixef.means[[k]]
  icm_k <- fit$fixef.mode[[k]]
  for (nm in names(dm_k)) {
    cat(sprintf("  %-28s  draws = %8.4f  ICM = %8.4f\n", nm, dm_k[[nm]], icm_k[[nm]]))
  }
}

## Demo: lmerb() with dGamma() observation dispersion on big_word_club
##
## Gaussian LMM route rLMMindepNormalGamma_reg_known_vcov(): random sigma^2
## (dGamma dispersion_ranef) with fixed Block~2 tau^2_k (all dNormal pfamily).
## Same model as demo/Ex_12_lmerb_BigWordClub.R; compare Ex_23 case 3
## (simulate = FALSE joint mode only).
##
## TEMP: subsets to algebraically full-rank schools. Also excludes school_id
## 2 and 18 (legacy per-group ING envelope issue). Block~1 uses the centering
## bridge (BlockEnvelopeCentering) until BlockEnvelopeSim ships; remove school
## filters once shared sigma^2 envelope is stable.
##
##   demo("Ex_24_lmerb_dGamma_BigWordClub", package = "lmebayes")

if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("This demo requires the 'bayesrules' package.", call. = FALSE)
}
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This demo requires the 'lme4' package.", call. = FALSE)
}

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

cat("\n=== lmer reference for BlockEnvelopeCentering check ===\n")
cat(sprintf(
  "  REML sigma^2 (compare to center$dispersion): %.6g\n",
  stats::sigma(fit_lmer)^2
))

grp_col      <- design$group_name
re_names_ref <- design$re_coef_names
fe_lmer      <- lme4::fixef(fit_lmer)
coef_raw     <- as.data.frame(coef(fit_lmer)[[grp_col]])
cn_lmer      <- names(coef_raw)
if (length(re_names_ref) == ncol(coef_raw)) {
  if (!is.null(cn_lmer) && !identical(cn_lmer, re_names_ref)) {
    if (all(re_names_ref %in% cn_lmer)) {
      coef_raw <- coef_raw[, re_names_ref, drop = FALSE]
    } else {
      names(coef_raw) <- re_names_ref
    }
  } else if (is.null(cn_lmer)) {
    names(coef_raw) <- re_names_ref
  }
}

cat("  X_hyper predictors per RE (mu_all = X_hyper %*% fixef from lmer fixef):\n")
for (k in re_names_ref) {
  cat(sprintf(
    "    %-20s: %s\n",
    k,
    paste(colnames(design$X_hyper[[k]]), collapse = ", ")
  ))
}

.fe_name_for_lmer <- function(k, col, fe) {
  if (k == "(Intercept)") {
    if (col %in% names(fe)) col else NA_character_
  } else if (col == "(Intercept)") {
    if (k %in% names(fe)) k else NA_character_
  } else {
    cand <- c(paste0(col, ":", k), paste0(k, ":", col))
    hit  <- cand[cand %in% names(fe)]
    if (length(hit)) hit[1L] else NA_character_
  }
}

fixef_lmer <- lapply(re_names_ref, function(k) {
  cols_k <- colnames(design$X_hyper[[k]])
  fe_nms <- vapply(
    cols_k, .fe_name_for_lmer, character(1L), k = k, fe = fe_lmer
  )
  miss <- is.na(fe_nms) | !fe_nms %in% names(fe_lmer)
  if (any(miss)) {
    stop(
      "lmer fixef missing term(s) for X_hyper[[", k, "]]: ",
      paste(cols_k[miss], collapse = ", "),
      call. = FALSE
    )
  }
  mu_k <- vapply(fe_nms, function(nm) unname(fe_lmer[nm]), numeric(1L))
  names(mu_k) <- cols_k
  mu_k
})
names(fixef_lmer) <- re_names_ref

coef_anchor <- vapply(re_names_ref, function(k) {
  if (k == "(Intercept)") {
    unname(fe_lmer["(Intercept)"])
  } else if (k %in% names(fe_lmer)) {
    unname(fe_lmer[k])
  } else {
    0
  }
}, numeric(1L))

mu_all_lmer <- build_mu_all(design, fixef_lmer)$mu_all

grp_levs  <- rownames(coef_raw)
lmer_full <- matrix(
  NA_real_,
  nrow = length(grp_levs),
  ncol = length(re_names_ref),
  dimnames = list(grp_levs, re_names_ref)
)
for (j in seq_len(nrow(coef_raw))) {
  lev <- grp_levs[j]
  for (k in re_names_ref) {
    lmer_full[lev, k] <- mu_all_lmer[k, lev] +
      (unname(coef_raw[[k]][j]) - coef_anchor[k])
  }
}

cat(
  "  lmer_full = mu_all(lmer fixef) + (coef - anchor); ",
  "compare to center$b_post_mean (centering uses ICM mu_all, not this table):\n"
)
print(round(lmer_full, 4))

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
stopifnot(is.matrix(fit$fixef.dispersion))
stopifnot(
  all(is.finite(fit$fixef.dispersion)), all(fit$fixef.dispersion > 0),
  all(apply(fit$fixef.dispersion, 2L, stats::sd) == 0)
)
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

cat("\n=== Random effects: lmer reference vs lmerb chain mean ===\n\n")
cat("  Pre-fit lmer_full (build_mu_all from lmer fixef) printed above;\n")
cat("  post-fit comparison uses fit$fixef.mu (same mu_all as the sampler).\n\n")
lmebayes:::print_mer_bayes_re_compare(fit)

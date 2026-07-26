# Manual smoke test: same as group_dgamma_smoke_test_allrank.R, but drops the
# three schools flagged by dGamma_list() asymmetric-window diagnostics
# (R_lo = lo_pct_BLUP/lo_pct_OLS < 0.25) on the full all-rank fixture:
#   school_id 6, 33, 41
#
# Model:
#   score_ppvt ~ 1 + distracted_ppvt + (1 + distracted_ppvt || school_id)
#
# Per-group dispersion_ranef (dGamma_list()) requires dispformula =
# ~<group_name> (here ~school_id); lmerb() rejects it under the ~1 default.
#
#   cd .../lmebayes
#   Rscript data-raw/group_dgamma_smoke_test_allrank_drop6_33_41.R
#
# Override chain count without editing the script:
#   set LMEBAYES_N_DGAMMA=3000
#   Rscript data-raw/group_dgamma_smoke_test_allrank_drop6_33_41.R

source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_DGAMMA <- as.integer(Sys.getenv("LMEBAYES_N_DGAMMA", unset = "3000"))
if (!is.finite(N_DGAMMA) || N_DGAMMA < 1L) {
  stop("LMEBAYES_N_DGAMMA must be a positive integer.", call. = FALSE)
}

DROP_SCHOOLS <- c("6", "33", "41")

fx <- .prepare_small5_all_full_rank_manual()
dat  <- fx$dat
form <- fx$form

drop <- intersect(DROP_SCHOOLS, levels(dat$school_id))
if (length(drop)) {
  cat(sprintf(
    "\n=== Dropping %d asymmetric-window outlier school(s): %s ===\n\n",
    length(drop),
    paste(drop, collapse = ", ")
  ))
  dat <- subset(dat, !as.character(school_id) %in% drop)
  dat$school_id <- droplevels(dat$school_id)
}

design <- model_setup(form, data = dat)
group_levels <- levels(design$group)
n_schools <- length(group_levels)
p_re <- length(design$re_coef_names)

stopifnot(identical(design$re_coef_names, c("(Intercept)", "distracted_ppvt")))
stopifnot(all(design$re_rank))
stopifnot(n_schools >= 1L)
stopifnot(!any(DROP_SCHOOLS %in% group_levels))

max_disp_perc <- 0.8

cat(sprintf(
  "\n=== All-rank group dGamma smoke test (drop 6/33/41): %d schools, n = %d chains ===\n\n",
  n_schools, N_DGAMMA
))

source("tests/manual/_reference_mer_compare.R")
ref_fits <- .print_reference_mer_compare(
  form        = form,
  dat         = dat,
  group_name  = "school_id",
  dispformula = ~school_id
)
fit_lmer <- ref_fits$lmer

## --- Prior setup + per-group dGamma_list() --------------------------------
ps <- Prior_Setup_lmebayes(
  form,
  data            = dat,
  pwt             = 0.01,
  pwt_measurement = 0.1,
  max_disp_perc   = max_disp_perc,
  dispformula     = ~school_id
)
pf <- pfamily_list(ps)
disp_pf_list <- dGamma_list(ps)

ing_grp <- ps$ing_prior_measurement_group
sigma2_hat_j <- vapply(ing_grp, `[[`, 0, "sigma2_hat")
names(sigma2_hat_j) <- group_levels

guard_df <- data.frame(
  school     = group_levels,
  n_j        = vapply(ing_grp, `[[`, 0, "n_j"),
  n_prior    = round(vapply(ing_grp, `[[`, 0, "n_prior"), 3),
  n_combined = round(vapply(ing_grp, `[[`, 0, "n_combined"), 3),
  sigma2_hat = round(sigma2_hat_j, 2),
  shape_ING  = round(vapply(ing_grp, `[[`, 0, "shape_ING"), 3),
  rate_gamma = round(vapply(ing_grp, `[[`, 0, "rate_gamma"), 1),
  pwt_group  = round(vapply(ing_grp, `[[`, 0, "pwt_group"), 4),
  disp_lower = round(vapply(disp_pf_list, function(pf) pf$prior_list$disp_lower, 0), 2),
  disp_upper = round(vapply(disp_pf_list, function(pf) pf$prior_list$disp_upper, 0), 2),
  stringsAsFactors = FALSE
)
cat("\n=== Per-group dGamma_list() measurement dispersion priors ===\n\n")
print(guard_df)

stopifnot(all(vapply(
  disp_pf_list,
  function(pf) pf$prior_list$disp_upper >= pf$prior_list$disp_lower,
  logical(1L)
)))

n_combined_j <- vapply(ing_grp, `[[`, 0, "n_combined")
rate_w_j <- sigma2_hat_j * (n_combined_j + p_re - 1) / 2
shape_w_j <- (n_combined_j + 1) / 2 + p_re / 2
disp_upper_sym <- 1 / stats::qgamma(1 - max_disp_perc, shape = shape_w_j, rate = rate_w_j)
stopifnot(all(vapply(disp_pf_list, function(pf) pf$prior_list$disp_upper, 0) >= disp_upper_sym - 1e-6))

cat("\n=== lmerb with per-group dGamma() priors; n =", N_DGAMMA, "===\n\n")

t_fit <- system.time({
  fit <- lmerb(
    form,
    data             = dat,
    pfamily_list     = pf,
    dispersion_ranef = disp_pf_list,
    dispformula      = ~school_id,
    n                = N_DGAMMA,
    progbar          = TRUE,
    verbose          = TRUE
  )
})
cat(sprintf("\n=== Timing: lmerb elapsed = %.2f s ===\n", t_fit["elapsed"]))

stopifnot(inherits(fit, "lmerb"))
stopifnot(identical(fit$prior$dispersion_mode, "gamma_list"))
stopifnot(identical(fit$draw_engine, "rGLMM_sweep_ing_block1_ind"))

stopifnot(is.matrix(fit$sigma2))
stopifnot(identical(dim(fit$sigma2), c(N_DGAMMA, n_schools)))
stopifnot(identical(colnames(fit$sigma2), group_levels))
stopifnot(length(fit$sigma2.mean) == n_schools)
stopifnot(identical(names(fit$sigma2.mean), group_levels))
stopifnot(all(is.finite(fit$sigma2)))
stopifnot(is.finite(fit$m_convergence) && fit$m_convergence >= 1L)

cat("\n=== Per-group sigma2 posterior means (should track sigma2_hat_j) ===\n\n")
print(data.frame(
  school       = group_levels,
  sigma2_hat_j = round(sigma2_hat_j[group_levels], 2),
  sigma2_mean  = round(fit$sigma2.mean[group_levels], 2)
))

cat(sprintf("\nm_convergence = %d, n_pilot = %s\n", fit$m_convergence, fit$convergence$n_pilot))

coef_focus <- list(
  c("(Intercept)", "(Intercept)"),
  c("distracted_ppvt", "(Intercept)")
)
cat("\n=== Block~2 convergence (plot_sweep_history_diag) ===\n\n")
if (!interactive()) {
  grDevices::pdf(NULL)
}
for (st in list(fit$sweep_history$pilot, fit$sweep_history$main)) {
  if (is.null(st)) next
  plot_sweep_history_diag(st, coef_focus)
}
if (!interactive() && grDevices::dev.cur() > 1L) {
  grDevices::dev.off()
}

cat("\n=== summary(lmerb) ===\n\n")
print(summary(fit))

cat("\n=== summary_sigma2(lmerb) — per-group measurement dispersion ===\n\n")
print(summary_sigma2(fit))

stopifnot(isTRUE(all.equal(lme4::fixef(fit$lmer), lme4::fixef(fit_lmer))))
stopifnot(isTRUE(all.equal(as.numeric(logLik(fit$lmer)), as.numeric(logLik(fit_lmer)))))

## --- error-path checks -----------------------------------------------------
if (n_schools >= 2L) {
  err <- tryCatch({
    lmerb(
      form, data = dat, pfamily_list = pf,
      dispersion_ranef = disp_pf_list[seq_len(n_schools - 1L)], n = 2L
    )
    NULL
  }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(err))
  cat("\n[expected error, wrong list length]:\n", err, "\n")
}

bad_list <- disp_pf_list
bad_list[[1L]]$prior_list$disp_lower <- NULL
bad_list[[1L]]$prior_list$disp_upper <- NULL
err2 <- tryCatch({
  lmerb(
    form, data = dat, pfamily_list = pf,
    dispersion_ranef = bad_list, n = 2L
  )
  NULL
}, error = function(e) conditionMessage(e))
stopifnot(!is.null(err2))
cat("\n[expected error, missing disp_lower/upper]:\n", err2, "\n")

cat("\ngroup_dgamma_smoke_test_allrank_drop6_33_41: OK\n")

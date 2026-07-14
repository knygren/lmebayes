# Manual smoke test: a NAMED LIST of dGamma() measurement-dispersion priors,
# one per group (the third `dispersion_ranef` option alongside a fixed
# scalar and a single pooled dGamma()).  Reuses the small5 fixture (5
# full-rank schools) so per-group draws can be sanity checked directly.
#
# Per-group dispersion_ranef (dGamma_list()) requires dispformula =
# ~<group_name> (here ~school_id); lmerb() rejects it under the ~1 default.
#
#   cd .../lmebayes
#   Rscript data-raw/group_dgamma_smoke_test.R

source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

N_DGAMMA <- 200L

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
dat  <- fx$dat
form <- fx$form
design <- fx$design
group_levels <- levels(design$groups)
p_re <- length(design$re_coef_names)

stopifnot(identical(design$re_coef_names, c("(Intercept)", "distracted_ppvt")))
stopifnot(all(design$re_rank))

max_disp_perc <- 0.99

## --- lme4 reference (tau^2 for diagnostics) --------------------------------
fit_lmer <- lme4::lmer(form, data = dat)

## --- Prior setup + per-group dGamma_list() --------------------------------
ps <- Prior_Setup_lmebayes(
  form,
  data            = dat,
  pwt             = 0.01,
  pwt_measurement = 0.1,
  max_disp_perc   = max_disp_perc
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

disp_pf_list

cat("\n=== lme4::lmer reference (same form and dat) ===\n\n")
print(summary(fit_lmer))
cat("\nfixef(lmer):\n")
print(lme4::fixef(fit_lmer))
cat("\nVarCorr(lmer):\n")
print(lme4::VarCorr(fit_lmer), comp = "Std.Dev.")
cat("\nranef(lmer):\n")
print(lme4::ranef(fit_lmer))
cat("\ncoef(lmer):\n")
print(coef(fit_lmer))

cat("\n=== lmerb with a per-group list of dGamma() priors; n =", N_DGAMMA, "===\n\n")

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
stopifnot(identical(dim(fit$sigma2), c(N_DGAMMA, 5L)))
stopifnot(identical(colnames(fit$sigma2), group_levels))
stopifnot(length(fit$sigma2.mean) == 5L)
stopifnot(identical(names(fit$sigma2.mean), group_levels))
stopifnot(all(is.finite(fit$sigma2)))
stopifnot(is.finite(fit$m_convergence) && fit$m_convergence >= 1L)

cat("\n=== Per-group sigma2 posterior means (should track sigma2_hat_j) ===\n\n")
print(data.frame(
  school        = group_levels,
  sigma2_hat_j  = round(sigma2_hat_j[group_levels], 2),
  sigma2_mean   = round(fit$sigma2.mean[group_levels], 2)
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

cat("\n=== lmerb extractors (lme4-style) ===\n\n")
cat("ranef(lmerb, type = 'mean', condVar = TRUE):\n")
print(lme4::ranef(fit, type = "mean", condVar = TRUE))
cat("\ncoef(lmerb, type = 'mean'):\n")
print(coef(fit, type = "mean"))
cat("\nfixef(lmerb) (= fixef(lmer)):\n")
print(lme4::fixef(fit))

stopifnot(isTRUE(all.equal(lme4::fixef(fit$lmer), lme4::fixef(fit_lmer))))
stopifnot(isTRUE(all.equal(as.numeric(logLik(fit$lmer)), as.numeric(logLik(fit_lmer)))))

## --- error-path checks -----------------------------------------------------
err <- tryCatch({
  lmerb(
    form, data = dat, pfamily_list = pf,
    dispersion_ranef = disp_pf_list[1:4], n = 2L
  )
  NULL
}, error = function(e) conditionMessage(e))
stopifnot(!is.null(err))
cat("\n[expected error, wrong list length]:\n", err, "\n")

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

cat("\ngroup_dgamma_smoke_test: OK\n")

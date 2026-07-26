# Isolate the per-group dGamma() measurement-dispersion prior (as used by the
# lmerb() gamma_list Block~1 path) from the confounding state of the full
# two-block Gibbs sweep, using lmb()/lmbBlock() -- which run the EXACT same
# rindepNormalGamma_reg() engine, one call per school, with a prior that
# looks AS CLOSE AS POSSIBLE to what the real Block~1 sampler actually uses:
#
#   - Sigma (prior covariance for b_j): diag(tau2), tau2 from ONE lmer() fit
#     on the pooled data (VarCorr) -- this IS exactly
#     glmbayesCore:::.lmebayes_block1_prior_list()'s `Sigma_ranef` (see
#     mixed_rmerb_helpers.R:370), not an sd-derived approximation with
#     spurious off-diagonal terms.
#   - mu (prior mean for b_j): fixef(lmer_fit) -- since X_hyper is
#     intercept-only here, build_mu_all() gives mu_all[, j] = fixef for
#     EVERY group j (see build_mu_all.R).
#   - shape/rate/disp_lower/disp_upper (dispersion prior): Prior_Setup()
#     per school with scalar pwt = 0.25 (already validated stable/well
#     accepted with a LOOSE coefficient prior in an earlier run of this
#     script; here we keep that dispersion calibration but swap in the
#     faithful mu/Sigma above).
#
#   cd .../lmebayes
#   Rscript data-raw/test_group_priors_via_lmbBlock.R

source("tests/manual/_load.R")
source("tests/manual/_small5_lmerb_fixture.R")
.manual_test_load(load_glmbayes_core = TRUE)

fx <- .prepare_small5_lmerb_manual(n_schools = 5L)
dat  <- fx$dat
design <- fx$design
group_levels <- levels(design$group)

## Per-block formula: same fixed-form as the RE design Z_j (Intercept + slope).
form_block <- score_ppvt ~ 1 + distracted_ppvt
form_lmer  <- score_ppvt ~ 1 + distracted_ppvt + (1 + distracted_ppvt || school_id)

n_j <- as.integer(table(design$group))
names(n_j) <- group_levels

## --- Sigma (diag(tau2)) and mu (fixef) from ONE lmer() fit -----------------
fit_lmer <- lme4::lmer(form_lmer, data = dat, REML = TRUE)
vc       <- as.data.frame(lme4::VarCorr(fit_lmer))
tau2 <- c(
  `(Intercept)`     = vc$vcov[vc$grp == "school_id"   & vc$var1 == "(Intercept)"],
  distracted_ppvt   = vc$vcov[vc$grp == "school_id.1" & vc$var1 == "distracted_ppvt"]
)
Sigma_re <- diag(tau2, nrow = 2L, ncol = 2L, names = TRUE)
dimnames(Sigma_re) <- list(names(tau2), names(tau2))
mu_re <- lme4::fixef(fit_lmer)[names(tau2)]

cat("\ntau2 (from lmer VarCorr):\n"); print(tau2)
cat("\nSigma_re = diag(tau2):\n"); print(Sigma_re)
cat("\nmu_re (from lmer fixef; same for every group):\n"); print(mu_re)

## --- Dispersion prior: Prior_Setup() per school, scalar pwt = 0.25 --------
pwt_measurement_group <- 0.25
ps_block <- stats::setNames(lapply(group_levels, function(lev) {
  dat_j <- dat[dat$school_id == lev, , drop = FALSE]
  Prior_Setup(
    form_block,
    data   = dat_j,
    family = gaussian(),
    pwt    = pwt_measurement_group
  )
}), group_levels)

## Smaller max_disp_perc => tighter posterior-connected dispersion bounds
## (see rationale below, at the dIndependent_Normal_Gamma() call).
max_disp_perc <- 0.95
guard_df <- data.frame(
  school     = group_levels,
  n_j        = n_j,
  sigma2_hat = round(vapply(ps_block, `[[`, 0, "dispersion"), 2),
  shape_ING  = round(vapply(ps_block, `[[`, 0, "shape_ING"), 3),
  rate_gamma = round(vapply(ps_block, `[[`, 0, "rate_gamma"), 1)
)
cat("\n=== Per-school dispersion prior (Prior_Setup(), pwt = 0.25) ===\n\n")
print(guard_df)

## --- Assemble dIndependent_Normal_Gamma() with FAITHFUL mu/Sigma and the
## validated dispersion shape/rate. disp_lower/disp_upper are left NULL: per
## EnvelopeDispersionBuild.cpp (Step 2), when NULL the C++ engine derives the
## dispersion bounds from the *posterior* Gamma(shape + n_w/2, rate +
## RSS_post/2) at max_disp_perc -- i.e. bounds connected to the posterior,
## not the (much wider) prior-only quantiles we were computing by hand
## before. Only max_disp_perc controls how tight those bounds are.
pfamily_list_j <- stats::setNames(lapply(group_levels, function(lev) {
  ps <- ps_block[[lev]]
  glmbayesCore::dIndependent_Normal_Gamma(
    mu            = mu_re,
    Sigma         = Sigma_re,
    shape         = ps$shape_ING,
    rate          = ps$rate_gamma,
    max_disp_perc = max_disp_perc
  )
}), group_levels)

N_TEST <- 200L

## --- Primary check: all 5 schools together via lmbBlock() -----------------
cat(sprintf("\n=== lmbBlock() across all %d schools (n = %d each) ===\n\n",
            length(group_levels), N_TEST))

out_blmb <- tryCatch(
  lmbBlock(
    form_block,
    block        = "school_id",
    pfamily_list = pfamily_list_j,
    data         = dat,
    n            = N_TEST,
    use_parallel = FALSE,
    verbose      = FALSE
  ),
  error = function(e) e
)

if (inherits(out_blmb, "error")) {
  cat("lmbBlock() FAILED:", conditionMessage(out_blmb), "\n")
  cat("Falling back to per-school isolation below.\n\n")
} else {
  print(out_blmb)
  cat("lmbBlock(): OK across all schools.\n\n")
}

## --- Per-school isolation (candidates/draw + explicit crash isolation) ----
cat("=== Per-school isolation (lmb() one school at a time) ===\n\n")
for (lev in group_levels) {
  dat_j <- dat[dat$school_id == lev, , drop = FALSE]
  cat(sprintf("--- school %s (n_j = %d) ---\n", lev, n_j[[lev]]))
  t0 <- system.time({
    fit_j <- tryCatch(
      lmb(
        form_block,
        pfamily      = pfamily_list_j[[lev]],
        data         = dat_j,
        n            = N_TEST,
        use_parallel = FALSE,
        verbose      = FALSE
      ),
      error = function(e) e
    )
  })
  if (inherits(fit_j, "error")) {
    cat("  FAILED:", conditionMessage(fit_j), "\n\n")
  } else {
    cat(sprintf(
      "  OK: mean(iters) = %.1f candidates/draw, mean(dispersion) = %.2f, elapsed = %.2fs\n\n",
      mean(fit_j$iters), mean(fit_j$dispersion), t0["elapsed"]
    ))
  }
}

cat("test_group_priors_via_lmbBlock: done\n")

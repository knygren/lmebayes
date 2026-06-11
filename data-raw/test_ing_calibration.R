# Regression test: dIndependent_Normal_Gamma pfamilies in lmerb()/glmerb().
#
# ING components require a positive 'disp_lower' (lower dispersion
# truncation), which is used as the conservative tau^2 plug-in for the
# eigenvalue / TV calibration: smaller tau^2 increases the block coupling
# and hence lambda*, so the disp_lower-based lambda* upper-bounds the rate
# for every dispersion in the truncated support.  Block 2 dispersion
# sampling is not implemented, so the fit displays the calibration and
# stops (no draws).
#
# Test value per the design discussion: disp_lower = tau^2_k / 2 (half the
# classical lmer-estimated RE variance).
#
#   Rscript data-raw/test_ing_calibration.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

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

ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
re_names <- names(ps$prior_list)

## --- baseline: dNormal path (full sampler) ----------------------------------
out_dn <- capture.output(
  fit_dn <- lmerb(form_lmer, data = dat,
                  pfamily_list = pfamily_list(ps),
                  dispersion_ranef = ps$dispersion_ranef,
                  n = 10L, seed = 1L)
)
stopifnot(
  !is.null(fit_dn$coefficients),
  identical(fit_dn$convergence$method, "exact")
)
lambda_dn <- fit_dn$convergence$lambda_star
m_min_dn  <- fit_dn$convergence$m_min

## --- ING pfamilies with disp_lower = tau^2_k / 2 ----------------------------
pf_ing0 <- pfamily_list(ps, ptypes = "dIndependent_Normal_Gamma")
pf_ing <- stats::setNames(lapply(re_names, function(k) {
  pr <- pf_ing0[[k]]$prior_list
  dIndependent_Normal_Gamma(
    mu         = pr$mu,
    Sigma      = pr$Sigma,
    shape      = pr$shape,
    rate       = pr$rate,
    disp_lower = unname(ps$prior_list[[k]]$dispersion_fixef) / 2
  )
}), re_names)

out_ing <- capture.output(
  fit_ing <- lmerb(form_lmer, data = dat,
                   pfamily_list = pf_ing,
                   dispersion_ranef = ps$dispersion_ranef,
                   n = 10L, seed = 1L)
)

## 1. Stops after calibration: no draws, but calibration info present/printed.
stopifnot(
  inherits(fit_ing, "lmerb"),
  is.null(fit_ing$coefficients),
  is.null(fit_ing$fixef_draws),
  is.null(fit_ing$coef.means),
  !is.null(fit_ing$coef.mode),
  !is.null(fit_ing$convergence),
  identical(fit_ing$convergence$method, "disp_lower_bound")
)
stopifnot(
  any(grepl("conservative: ING tau\\^2_k = disp_lower", out_ing)),
  any(grepl("Stopping after the", out_ing))
)
cat(sprintf(
  "1. lmerb ING: stopped after calibration; lambda* = %.4f (dNormal: %.4f), m_min = %d (dNormal: %d)\n",
  fit_ing$convergence$lambda_star, lambda_dn,
  fit_ing$convergence$m_min, m_min_dn
))

## 2. Conservative ordering: halving tau^2 increases coupling => larger
##    lambda* and at least as many sweeps.
stopifnot(
  fit_ing$convergence$lambda_star > lambda_dn,
  fit_ing$convergence$m_min >= m_min_dn
)
cat("2. lambda*(disp_lower = tau^2/2) > lambda*(tau^2): OK\n")

## 3. The plug-in matches an explicit dNormal fit at tau^2/2 exactly.
pf_half <- stats::setNames(lapply(re_names, function(k) {
  pl <- ps$prior_list[[k]]
  dNormal(mu = pl$mu_fixef, Sigma = pl$Sigma_fixef,
          dispersion = unname(pl$dispersion_fixef) / 2)
}), re_names)
out_half <- capture.output(
  fit_half <- lmerb(form_lmer, data = dat,
                    pfamily_list = pf_half,
                    dispersion_ranef = ps$dispersion_ranef,
                    n = 10L, seed = 1L)
)
stopifnot(isTRUE(all.equal(
  fit_ing$convergence$lambda_star, fit_half$convergence$lambda_star
)))
stopifnot(identical(fit_ing$convergence$m_min, fit_half$convergence$m_min))
cat("3. ING calibration == dNormal calibration at tau^2/2: OK\n")

## 4. Missing disp_lower errors clearly (constructed without one; the
##    pfamily_list() builder now sets a default, so build manually).
pf_bad <- pf_ing
pr1 <- pf_ing0[[1L]]$prior_list
pf_bad[[1L]] <- dIndependent_Normal_Gamma(
  mu = pr1$mu, Sigma = pr1$Sigma, shape = pr1$shape, rate = pr1$rate
)
err <- tryCatch(
  lmerb(form_lmer, data = dat, pfamily_list = pf_bad,
        dispersion_ranef = ps$dispersion_ranef, n = 5L),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err), grepl("disp_lower", err))
cat("4. Missing disp_lower rejected: OK\n")

## 4b. pfamily_list() default disp_lower (0.01 dispersion quantile =
##     1/qgamma(0.99, shape, rate)) passes lmerb validation end-to-end and
##     stops after calibration.  With the diffuse default calibration the
##     quantile sits far below tau^2, so lambda* should exceed the tau^2/2
##     case.
for (k in re_names) {
  pr <- pf_ing0[[k]]$prior_list
  stopifnot(isTRUE(all.equal(
    pr$disp_lower, 1 / stats::qgamma(0.99, shape = pr$shape, rate = pr$rate)
  )))
}
out_def <- capture.output(
  fit_def <- lmerb(form_lmer, data = dat,
                   pfamily_list = pf_ing0,
                   dispersion_ranef = ps$dispersion_ranef,
                   n = 10L, seed = 1L)
)
stopifnot(
  is.null(fit_def$coefficients),
  identical(fit_def$convergence$method, "disp_lower_bound"),
  fit_def$convergence$lambda_star > fit_ing$convergence$lambda_star,
  fit_def$convergence$m_min >= fit_ing$convergence$m_min
)
ratios <- vapply(re_names, function(k) {
  pf_ing0[[k]]$prior_list$disp_lower /
    unname(ps$prior_list[[k]]$dispersion_fixef)
}, numeric(1L))
cat(sprintf(
  "4b. default disp_lower: lambda* = %.4f, m_min = %d (disp_lower/tau^2 ratios: %s)\n",
  fit_def$convergence$lambda_star, fit_def$convergence$m_min,
  paste(sprintf("%.3f", ratios), collapse = ", ")
))

## 5. Mixed list (ING + dNormal) also stops after calibration.
pf_mixed <- pfamily_list(ps)
pf_mixed[[2L]] <- pf_ing[[2L]]
out_mix <- capture.output(
  fit_mix <- lmerb(form_lmer, data = dat,
                   pfamily_list = pf_mixed,
                   dispersion_ranef = ps$dispersion_ranef,
                   n = 10L, seed = 1L)
)
stopifnot(
  is.null(fit_mix$coefficients),
  identical(fit_mix$convergence$method, "disp_lower_bound"),
  fit_mix$convergence$lambda_star > lambda_dn
)
cat("5. Mixed dNormal/ING list stops after calibration: OK\n")

## 6. glmerb (gaussian) path: same behavior, combined method label.
out_g <- capture.output(
  fit_g <- glmerb(form_lmer, data = dat, family = gaussian(),
                  pfamily_list = pf_ing,
                  dispersion_ranef = ps$dispersion_ranef,
                  n = 10L, seed = 1L)
)
stopifnot(
  inherits(fit_g, "glmerb"),
  is.null(fit_g$coefficients),
  identical(fit_g$convergence$method, "exact+disp_lower_bound"),
  isTRUE(all.equal(fit_g$convergence$lambda_star,
                   fit_ing$convergence$lambda_star))
)
stopifnot(any(grepl("conservative: ING tau\\^2_k = disp_lower", out_g)))
cat("6. glmerb gaussian ING path: OK\n")

cat("\ntest_ing_calibration.R: all checks passed\n")

## Prior_Setup_lmebayes -- development script
##
## Drafts the Prior_Setup_lmebayes() function for the two-block Gibbs sampler.
##
## First implementation: Gaussian family, variance components FIXED at lmer
## estimates (dNormal analog from glmbayes).  The posterior over the unknowns
## is therefore a single large multivariate normal, sampled in two blocks:
##
##   Block 1  (per-group, independent across j):
##     p(b[j] | y, fixef, dispersion_ranef, Sigma_ranef)
##       = N(mu_b_post[j], Sigma_b_post[j])
##
##       Sigma_b_post[j]^{-1} = Z_j'Z_j / dispersion_ranef
##                              + diag(1 / diag(Sigma_ranef))
##       mu_b_post[j]          = Sigma_b_post[j] *
##                               (Z_j'y[j] / dispersion_ranef
##                                + diag(1/diag(Sigma_ranef)) %*% mu_b_prior[j])
##
##     where mu_b_prior[j][k] = X_hyper_k[j,] %*% prior_list[[k]]$mu_fixef
##
##   Block 2  (per-RE k, independent across k):
##     p(fixef_k | b_k, dispersion_fixef_k)
##       = N(mu_fixef_post_k, Sigma_fixef_post_k)
##
##       Sigma_fixef_post_k^{-1} = X_k'X_k / dispersion_fixef_k
##                                 + Sigma_fixef_k^{-1}
##       mu_fixef_post_k          = Sigma_fixef_post_k *
##                                  (X_k'b_k / dispersion_fixef_k
##                                   + Sigma_fixef_k^{-1} %*% mu_fixef_k)
##
##     where b_k = (b_k[1], ..., b_k[J])  current per-group slopes for RE k
##
## Prior returned by Prior_Setup_lmebayes():
##   dispersion_ranef : sigma2  (scalar, fixed at lmer estimate)
##   Sigma_ranef      : diagonal matrix (p_re x p_re), tau2_k on diagonal
##   prior_list       : named list, one entry per RE k, each containing:
##                        $mu_fixef        -- prior mean vector for fixef_k
##                        $Sigma_fixef     -- prior covariance matrix for fixef_k
##                        $dispersion_fixef-- tau2_k scalar (= Sigma_ranef[k,k])
##
## Calibration philosophy:
##   mu_fixef_k    <- fixef(fit_fr) for every X_hyper column (each random slope
##                    needs a matching fixed main effect in the formula)
##   Sigma_fixef_k <- vcov(fit_fr)[relevant cols] * (1-pwt)/pwt
##   tau2_k        <- VarCorr(fit_fr); must be strictly positive (non-singular)
##   Scaling (1-pwt)/pwt matches glmbayes::compute_gaussian_prior which uses
##   Sigma = (n_eff/n_prior) * dispersion * (X'X)^{-1} and n_eff/n_prior = (1-pwt)/pwt.
##
## All quantities are calibrated from the full-rank-group refit so that
## rank-deficient groups (zero-imputed BLUPs) do not bias the prior.

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4: install.packages('lme4')")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules: install.packages('bayesrules')")
}
pkgload::load_all(export_all = FALSE)

## Prior_Setup_lmebayes() and print.lmebayes_prior_setup are now in
## R/prior_setup_lmebayes.R.  This script loads them via pkgload::load_all().

## ===========================================================================
## Development run: big_word_club (same formula as inst/examples/Ex_lmerb.R)
## ===========================================================================

data(big_word_club, package = "bayesrules")
dat <- big_word_club
dat$school_id <- factor(dat$school_id)
dat <- subset(
  dat,
  !is.na(score_ppvt) &
    !is.na(invalid_ppvt) & invalid_ppvt == 0L &
    complete.cases(dat[, c("score_ppvt", "distracted_a1", "distracted_ppvt",
                           "private_school", "title1", "free_reduced_lunch",
                           "school_id")])
)

form_lmer <- score_ppvt ~
  private_school + title1 + free_reduced_lunch +
  distracted_a1 + distracted_ppvt +
  free_reduced_lunch:distracted_a1 +
  (1 + distracted_ppvt + distracted_a1 || school_id)

cat("=== Prior_Setup_lmebayes: pwt = 0.01 ===\n\n")
ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
print(ps)

## ---------------------------------------------------------------------------
## Sanity check: mu_fixef matches clean lmer fixed-effect estimates
## ---------------------------------------------------------------------------
cat("\n=== Sanity check: mu_fixef matches lmer fixed-effect estimates ===\n\n")

for (nm in ps$design$re_coef_names) {
  pl   <- ps$prior_list[[nm]]
  cols <- names(pl$mu_fixef)
  cat(sprintf("  [%s]  (hyper model: ~ %s)\n", nm,
              if (length(cols) == 1L) "1"
              else paste(c("1", cols[-1L]), collapse = " + ")))
  cat(sprintf("    %-22s  %s\n", "X_hyper column", "mu_fixef"))
  for (col in cols) {
    cat(sprintf("    %-22s  %10.6f\n", col, pl$mu_fixef[col]))
  }
  cat("\n")
}

## ---------------------------------------------------------------------------
## Block 1 demo: posterior b[j] for the first group
##   mu_b_prior[j][k] = X_hyper_k[j,] %*% mu_fixef_k
##   Sigma_b_post^{-1} = Z_j'Z_j / dispersion_ranef + Sigma_ranef^{-1}
##   mu_b_post         = Sigma_b_post * (Z_j'y_j / dispersion_ranef
##                                       + Sigma_ranef^{-1} %*% mu_b_prior)
## ---------------------------------------------------------------------------
cat("=== Block 1 demo: posterior b[j] for first group ===\n")

design   <- ps$design
grp1     <- levels(design$groups)[1L]
idx1     <- which(as.character(design$groups) == grp1)
y1       <- design$y[idx1]
Z1       <- design$Z[idx1, , drop = FALSE]
Sig_r    <- ps$Sigma_ranef
Sig_r_inv <- diag(1 / diag(Sig_r), nrow = nrow(Sig_r))   # diagonal inverse

mu_b_prior1 <- setNames(
  vapply(design$re_coef_names, function(k) {
    xrow <- design$X_hyper[[k]][grp1, , drop = TRUE]
    sum(xrow * ps$prior_list[[k]]$mu_fixef)
  }, numeric(1L)),
  design$re_coef_names
)

Sigma_b_post <- solve(crossprod(Z1) / ps$dispersion_ranef + Sig_r_inv)
mu_b_post    <- as.numeric(
  Sigma_b_post %*% (crossprod(Z1, y1) / ps$dispersion_ranef +
                      Sig_r_inv %*% mu_b_prior1)
)
names(mu_b_post) <- design$re_coef_names
dimnames(Sigma_b_post) <- list(design$re_coef_names, design$re_coef_names)

cat(sprintf("  group: %s   n_obs: %d\n\n", grp1, length(y1)))
cat("  mu_b_prior (X_hyper[j,] %*% mu_fixef):\n"); print(round(mu_b_prior1, 4L))
cat("\n  mu_b_post:\n");                            print(round(mu_b_post,   4L))
cat("\n  Sigma_b_post:\n");                         print(round(Sigma_b_post, 4L))
cat("\n  lmer ranef (BLUP, for reference):\n")
blup1 <- as.numeric(lme4::ranef(design$lmer_fit)[[design$group_name]][grp1, ])
print(round(setNames(blup1, design$re_coef_names), 4L))

## ---------------------------------------------------------------------------
## build_mu_all: full mu_all matrix for Block 1 (iter-0 fixef from prior_list)
## ---------------------------------------------------------------------------
cat("\n=== build_mu_all (iter 0) ===\n")
fixef0 <- lapply(ps$prior_list, `[[`, "mu_fixef")
names(fixef0) <- ps$design$re_coef_names
mu_out <- build_mu_all(ps$design, fixef0)
cat(sprintf("  mu_all: %d x %d (RE x groups)\n", nrow(mu_out$mu_all), ncol(mu_out$mu_all)))
cat("  first group, all RE:\n")
print(round(mu_out$mu_all[, 1L, drop = TRUE], 4L))
stopifnot(identical(mu_out$re_coef_names, ps$design$re_coef_names))
stopifnot(all.equal(
  as.numeric(mu_out$mu_all[, 1L]),
  as.numeric(mu_b_prior1),
  tolerance = 1e-8
))

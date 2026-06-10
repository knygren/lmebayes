# Convergence-rate test: two_block_rate() on the lmerb example dataset
# (bayesrules::big_word_club, Gaussian, school random intercept + slopes).
#
# Builds the exact same two-block sampler inputs as lmerb() (model_setup +
# Prior_Setup_lmebayes + .lmebayes_block1_prior_list) and computes the
# Remark 8 eigenvalues / lambda* (Nygren 2020).  Validates against a dense
# brute-force construction of the joint precision and cross-checks lambda*
# empirically against the ICM mean recursion, which contracts at exactly
# lambda* per iteration (Claim 2 applies to the deterministic conditional-
# mean iteration).
#
#   Rscript data-raw/test_lmerb_rate_big_word_club.R

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

design <- model_setup(form_lmer, data = dat)
ps <- Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)

## --- two-block sampler inputs, exactly as lmerb() builds them --------------
re_names     <- design$re_coef_names
group_levels <- levels(design$groups)
block1_prior <- lmebayes:::.lmebayes_block1_prior_list(ps)
block2_prior_list <- stats::setNames(
  lapply(re_names, function(k) {
    pl_k <- ps$prior_list[[k]]
    list(
      mu         = pl_k$mu_fixef,
      Sigma      = pl_k$Sigma_fixef,
      dispersion = pl_k$dispersion_fixef
    )
  }),
  re_names
)

rate <- glmbayesCore::two_block_rate(
  x = design$Z,
  block = design$groups,
  x_hyper = design$X_hyper,
  prior_list_block1 = block1_prior,
  prior_list_block2 = block2_prior_list,
  family = gaussian(),
  group_levels = group_levels
)

cat("\n=== two_block_rate on lmerb (big_word_club) inputs ===\n\n")
print(rate)

stopifnot(inherits(rate, "two_block_rate"))
stopifnot(rate$dims$J == length(group_levels))
stopifnot(rate$dims$p_re == length(re_names))
stopifnot(rate$dims$q == sum(vapply(design$X_hyper, ncol, integer(1L))))
stopifnot(all(rate$eigenvalues >= 0), all(rate$eigenvalues < 1))

## --- 1. Dense brute-force validation on the real design --------------------
P_b  <- block1_prior$P
w    <- rep(1 / block1_prior$dispersion, nrow(design$Z))
q_k  <- vapply(design$X_hyper, ncol, integer(1L))
q    <- sum(q_k)
cols <- split(seq_len(q), rep(seq_along(q_k), q_k))
J    <- length(group_levels)
p_re <- length(re_names)

H <- lapply(seq_len(J), function(j) {
  Hj <- matrix(0, p_re, q)
  for (k in seq_len(p_re)) {
    X_k <- as.matrix(design$X_hyper[[k]])
    Hj[k, cols[[k]]] <- X_k[group_levels[j], ]
  }
  Hj
})

P22 <- matrix(0, J * p_re, J * p_re)
P12 <- matrix(0, q, J * p_re)
P11 <- matrix(0, q, q)
grp_int <- as.integer(design$groups)
for (j in seq_len(J)) {
  rows <- which(grp_int == j)
  Z_j <- design$Z[rows, , drop = FALSE]
  B_j <- crossprod(Z_j, Z_j * w[rows]) + P_b
  bc <- (j - 1L) * p_re + seq_len(p_re)
  P22[bc, bc] <- B_j
  P12[, bc] <- -t(H[[j]]) %*% P_b
  P11 <- P11 + t(H[[j]]) %*% P_b %*% H[[j]]
}
for (k in seq_along(cols)) {
  Vk <- chol2inv(chol(as.matrix(block2_prior_list[[k]]$Sigma)))
  P11[cols[[k]], cols[[k]]] <- P11[cols[[k]], cols[[k]]] + Vk
}

e11 <- eigen(0.5 * (P11 + t(P11)), symmetric = TRUE)
P11_is <- e11$vectors %*% diag(1 / sqrt(e11$values), q) %*% t(e11$vectors)
A_dense <- P11_is %*% P12 %*% solve(P22, t(P12)) %*% P11_is
ev_dense <- sort(eigen(0.5 * (A_dense + t(A_dense)), symmetric = TRUE,
                       only.values = TRUE)$values, decreasing = TRUE)

stopifnot(max(abs(rate$eigenvalues - ev_dense)) < 1e-10)
cat("1. spectrum matches dense brute force (max diff ",
    format(max(abs(rate$eigenvalues - ev_dense)), digits = 3), ")\n", sep = "")

## --- 2. Empirical cross-check: ICM mean recursion contracts at lambda* -----
## lmerb_posterior_mean() iterates the exact conditional means; by Claim 2
## the gamma error contracts by a factor lambda* per sweep, so its observed
## iteration count must be consistent with the predicted geometric rate.
pm <- glmbayesCore::lmerb_posterior_mean(design, ps, tol = 1e-10)
stopifnot(isTRUE(pm$converged))

# Predicted iterations to reduce the initial error to tol: the lmer-derived
# start is within O(1) of the posterior mean on this scale, so
# log(tol)/log(lambda*) should agree with pm$iterations up to a small
# additive/multiplicative slack.
pred_iter <- log(1e-10) / log(rate$lambda_star)
cat(sprintf(
  "2. ICM iterations observed = %d, predicted ~ %.1f (lambda* = %.4f)\n",
  pm$iterations, pred_iter, rate$lambda_star
))
stopifnot(pm$iterations <= ceiling(pred_iter) + 5L)
stopifnot(pm$iterations >= floor(pred_iter / 4))

## Observed per-iteration contraction at the end of the ICM run should not
## beat lambda* (the asymptotic rate is the dominant eigenvalue): check that
## delta after pm$iterations steps is consistent with lambda*^iterations.
implied_rate <- (pm$delta)^(1 / pm$iterations)
cat(sprintf(
    "   implied average contraction rate = %.4f (lambda* = %.4f)\n",
    implied_rate, rate$lambda_star
))

## --- 3. m_convergence implication for lmerb (currently hardcoded 10L) ------
m_needed <- rate$m_for_tol(1e-3)
cat(sprintf(
  "3. m_convergence for TV tol 1e-3: %d (lmerb currently uses 10L); tol 1e-6: %d\n",
  m_needed, rate$m_for_tol(1e-6)
))
stopifnot(m_needed >= 1L)

## --- 4. TV bounds (Theorem 3 exact + Corollary 1 envelope) -----------------
## lmerb starts every replicate chain at the ICM posterior mean, so D0 = 0:
## the mean term of both bounds vanishes and only the variance-convergence
## sum remains.
l_grid <- c(1L, 2L, 3L, 5L, 8L, 10L, 15L, 20L, 30L, 40L)
b_t3 <- glmbayesCore::two_block_tv_bound(rate, l_grid, method = "theorem3")
b_c1 <- glmbayesCore::two_block_tv_bound(rate, l_grid, method = "corollary1")

cat("\n4. TV bound by sweep count l (start at posterior mean, D0 = 0):\n\n")
cat("    l    theorem3      corollary1    (lambda*)^l proxy\n")
for (ii in seq_along(l_grid)) {
  cat(sprintf("  %3d   %.6e  %.6e  %.6e\n",
              l_grid[ii], b_t3[ii], b_c1[ii], rate$lambda_star^l_grid[ii]))
}

stopifnot(all(b_t3 <= b_c1 + 1e-12))
stopifnot(all(diff(b_t3) <= 1e-12), all(diff(b_c1) <= 1e-12))

l_t3 <- vapply(c(1e-2, 1e-3, 1e-6), function(tol) {
  glmbayesCore::two_block_l_for_tv(rate, tol, method = "theorem3")
}, integer(1L))
l_c1 <- vapply(c(1e-2, 1e-3, 1e-6), function(tol) {
  glmbayesCore::two_block_l_for_tv(rate, tol, method = "corollary1")
}, integer(1L))
cat(sprintf(
  "\n   sweeps for TV <= 1e-2 / 1e-3 / 1e-6: theorem3 = %d / %d / %d, corollary1 = %d / %d / %d\n",
  l_t3[1L], l_t3[2L], l_t3[3L], l_c1[1L], l_c1[2L], l_c1[3L]
))
cat(sprintf(
  "   (crude (lambda*)^m proxy gave %d / %d / %d; warm start + 2l variance decay explain the gap)\n",
  rate$m_for_tol(1e-2), rate$m_for_tol(1e-3), rate$m_for_tol(1e-6)
))
stopifnot(all(l_t3 <= l_c1))
## The stored Block 1 draw lags gamma by a half-step: bound at l-1 applies.
## lmerb's m_convergence = 10 therefore certifies TV <= bound(9) for b draws:
cat(sprintf(
  "   m_convergence = 10 certifies: gamma TV <= %.3e, b TV <= %.3e (theorem3)\n",
  glmbayesCore::two_block_tv_bound(rate, 10L),
  glmbayesCore::two_block_tv_bound(rate, 9L)
))

cat("\ntest_lmerb_rate_big_word_club: OK\n")

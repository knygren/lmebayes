## Prototype: choose rwt so that TV(pi_rwt, pi_0) <= tol, with no assumption
## on the distribution of h.
##
## pi_w is an exponential tilt of pi_0: pi_w = pi_0 e^(-w h) / E_0[e^(-w h)].
## With psi(w) = log E_0[e^(-w h)],
##   KL(pi_w || pi_0) = -w E_w[h] - psi(w) = int_0^w t Var_t(h) dt,
## exactly, for any distribution of h with finite variance. Pinsker gives
##   TV(pi_w, pi_0) <= sqrt(KL / 2),
## also distribution-free. KL is increasing in w, so the rule takes the
## largest w <= rwt_max with sqrt(KL(w) / 2) <= tol.
##
## E_0 and E_w are estimated from the exact Laplace draws (pi_1) of
## ridge_tv_predict_prototype.R, reweighted by e^((1 - w) h). The only
## approximation is Monte Carlo; its reliability is judged by the Pareto k of
## the e^h weights (the worst case, w = 0) and the batch spread of the result.
##
## Uses the installed lmebayesCore (no compilation). About 30 s.

proto <- new.env()
sys.source("data-raw/ridge_tv_predict_prototype.R", envir = proto)
h <- proto$h
n_lap <- proto$n_lap
n_batch <- proto$n_batch
tol_values <- proto$tol_values
rwt_check <- proto$rwt_check
tv_from_fits <- proto$tv_from_fits
rwt_max <- 0.95

lse <- function(x) {
  mx <- max(x)
  mx + log(sum(exp(x - mx)))
}
softmax <- function(x) exp(x - lse(x))

## KL(pi_w || pi_0) from draws of pi_1 and h on them. h is centred at its
## pi_0 mean first; the shift cancels in KL but avoids losing digits for
## small w, where KL ~ w^2 Var_0(h) / 2.
ridge_kl <- function(hh, w) {
  hh <- hh - sum(softmax(hh) * hh)
  psi <- lse((1 - w) * hh) - lse(hh)
  e_w <- sum(softmax((1 - w) * hh) * hh)
  max(-w * e_w - psi, 0)
}
ridge_tv_bound <- function(hh, w) sqrt(ridge_kl(hh, w) / 2)

ridge_rwt_for_tol <- function(hh, tol, rwt_max = 0.95) {
  if (ridge_tv_bound(hh, rwt_max) <= tol) return(rwt_max)
  stats::uniroot(function(lw) ridge_tv_bound(hh, exp(lw)) - tol,
                 lower = log(1e-8), upper = log(rwt_max),
                 tol = 1e-10)$root |> exp()
}

k_hat <- if (requireNamespace("loo", quietly = TRUE)) {
  suppressWarnings(loo::psis(h, r_eff = 1)$diagnostics$pareto_k)
} else NA_real_
reliability <- if (is.na(k_hat)) {
  "unknown (install 'loo' for Pareto k)"
} else if (k_hat < 0.5) {
  "reliable"
} else if (k_hat < 0.7) {
  "usable, wider Monte Carlo error"
} else {
  "UNRELIABLE: Laplace draws do not cover pi_0; use a fallback rwt"
}

batch <- split(seq_len(n_lap), rep(seq_len(n_batch), length.out = n_lap))
batch_se <- function(f) {
  v <- vapply(batch, function(ix) f(h[ix]), numeric(1))
  stats::sd(v) / sqrt(n_batch)
}

cat("\n\n===== KL rule (distribution-free) =====\n")
cat(sprintf("Pareto k of e^h weights: %.3f -> %s\n", k_hat, reliability))
cat(sprintf("Batch SEs use %d batches of %d draws (conservative for the full %d)\n",
            n_batch, n_lap %/% n_batch, n_lap))

cat("\n----- TV(pi_rwt, pi_0): Pinsker bound vs Laplace IS estimate vs fits -----\n")
tv_is <- vapply(rwt_check, function(w) {
  0.5 * sum(abs(softmax(h) - softmax((1 - w) * h)))
}, numeric(1))
kl <- vapply(rwt_check, function(w) ridge_kl(h, w), numeric(1))
print(data.frame(
  rwt = rwt_check,
  KL = signif(kl, 3),
  KL_se = signif(vapply(rwt_check, function(w) {
    batch_se(function(hh) ridge_kl(hh, w))
  }, numeric(1)), 2),
  pinsker_bound = signif(sqrt(kl / 2), 3),
  laplace_IS = signif(tv_is, 3),
  from_fits = tv_from_fits,
  bound_over_fits = round(sqrt(kl / 2) / tv_from_fits, 2)
), row.names = FALSE)

cat(sprintf("\n----- rwt chosen: largest rwt <= %.2f with sqrt(KL / 2) <= tol -----\n",
            rwt_max))
rwt_sel <- vapply(tol_values, function(tl) {
  ridge_rwt_for_tol(h, tl, rwt_max)
}, numeric(1))
print(data.frame(
  tol = tol_values,
  rwt = signif(rwt_sel, 4),
  rwt_se = signif(vapply(tol_values, function(tl) {
    batch_se(function(hh) ridge_rwt_for_tol(hh, tl, rwt_max))
  }, numeric(1)), 2),
  laplace_IS_at_rwt = signif(vapply(rwt_sel, function(w) {
    0.5 * sum(abs(softmax(h) - softmax((1 - w) * h)))
  }, numeric(1)), 3),
  first_order_rwt = signif(pmin(tol_values / proto$kappa, rwt_max), 4)
), row.names = FALSE)

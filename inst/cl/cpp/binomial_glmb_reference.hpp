// Reference helpers mirroring glmbayes::famfuncs_binomial.cpp dbinom_glmb():
// y_prop is the GLM binomial success proportion; wt is trial count (prior weight).
//
// Not part of the R package compile — for parity checks vs OpenCL kernels.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace glmbayes::cl_ref {

inline void binomial_glmb_round_counts(double y_prop, double wt, int* trials_out,
                                       int* success_out) {
  *trials_out = static_cast<int>(std::round(wt));
  *success_out = static_cast<int>(std::round(y_prop * wt));
}

inline double clamp_prob(double p) {
  return std::min(1.0, std::max(0.0, p));
}

/// Integer binomial log pmf; matches R::dbinom(k, n, p, TRUE) for valid integers.
inline double dbinom_int_logpmf(int k, int n, double p) {
  p = clamp_prob(p);
  const double q = 1.0 - p;
  if (n < 0) return -std::numeric_limits<double>::infinity();
  if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
  if (p <= 0.0) return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
  if (p >= 1.0) return (k == n) ? 0.0 : -std::numeric_limits<double>::infinity();
  const double lc =
      std::lgamma(static_cast<double>(n + 1)) -
      std::lgamma(static_cast<double>(k + 1)) -
      std::lgamma(static_cast<double>(n - k + 1));
  return lc + static_cast<double>(k) * std::log(p) +
         static_cast<double>(n - k) * std::log(q);
}

/// Per-row contribution matching yy = -dbinom_glmb(..., lg=TRUE) on CPU path.
inline double neg_ll_binomial_glmb(double y_prop, double wt, double mean_p_unclamped) {
  int ntr = 0, ks = 0;
  binomial_glmb_round_counts(y_prop, wt, &ntr, &ks);
  const double pc = clamp_prob(mean_p_unclamped);
  return -dbinom_int_logpmf(ks, ntr, pc);
}

}  // namespace glmbayes::cl_ref

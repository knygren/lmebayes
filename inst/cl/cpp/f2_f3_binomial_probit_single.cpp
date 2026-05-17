// Single grid-point reference for f2_f3_binomial_probit (matches inst/cl/src kernel layout).

#include "binomial_glmb_reference.hpp"
#include <cmath>
#include <cstring>

namespace glmbayes::cl_ref {

/// Standard normal cdf/pdf — reference only (kernel uses pnorm5/dnorm4).
inline double phi_cdf(double x) {
  return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

inline double phi_pdf(double x) {
  constexpr double kTwoPi = 6.283185307179586476925286766559;
  return std::exp(-0.5 * x * x) / std::sqrt(kTwoPi);
}

inline void f2_f3_binomial_probit_single(const double* X, const double* B,
                                         const double* mu, const double* P,
                                         const double* alpha, const double* y,
                                         const double* wt, int l1, int l2,
                                         double* qf_out, double* grad_out) {
  double tmp[64], g_loc[64];

  for (int k = 0; k < l2; ++k) {
    double acc = 0.0;
    for (int ell = 0; ell < l2; ++ell) {
      acc += P[k * l2 + ell] * (B[ell] - mu[ell]);
    }
    tmp[k] = acc;
  }

  double qsum = 0.0;
  for (int k = 0; k < l2; ++k) {
    qsum += (B[k] - mu[k]) * tmp[k];
  }
  double res_acc = 0.5 * qsum;
  std::memcpy(g_loc, tmp, sizeof(double) * static_cast<size_t>(l2));

  for (int i = 0; i < l1; ++i) {
    double eta = alpha[i];
    for (int k = 0; k < l2; ++k) {
      eta += X[k * l1 + i] * B[k];
    }

    const double p1 = phi_cdf(eta);
    const double p2 = phi_cdf(-eta);
    const double d = phi_pdf(eta);

    res_acc += neg_ll_binomial_glmb(y[i], wt[i], p1);

    const double resid = ((y[i] * d / p1) - ((1.0 - y[i]) * d / p2)) * wt[i];
    for (int k = 0; k < l2; ++k) {
      g_loc[k] -= X[k * l1 + i] * resid;
    }
  }

  *qf_out = res_acc;
  std::memcpy(grad_out, g_loc, sizeof(double) * static_cast<size_t>(l2));
}

}  // namespace glmbayes::cl_ref

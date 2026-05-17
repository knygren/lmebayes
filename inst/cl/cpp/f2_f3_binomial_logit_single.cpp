// Single grid-point reference for f2_f3_binomial_logit (matches inst/cl/src kernel layout).

#include "binomial_glmb_reference.hpp"
#include <cstring>

namespace glmbayes::cl_ref {

/// One OpenCL work-item: X column-major (index k*l1+i), B row-major coefficients,
/// P stored like kernel (index k*l2+ell), symmetric precision OK.
inline void f2_f3_binomial_logit_single(const double* X, const double* B,
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
    double dot = -alpha[i];
    for (int k = 0; k < l2; ++k) {
      dot -= X[k * l1 + i] * B[k];
    }
    double e, p, q;
    if (dot <= 0.0) {
      e = std::exp(dot);
      p = 1.0 / (1.0 + e);
      q = e / (1.0 + e);
    } else {
      e = std::exp(-dot);
      p = e / (1.0 + e);
      q = 1.0 / (1.0 + e);
    }

    res_acc += neg_ll_binomial_glmb(y[i], wt[i], p);

    const double resid = (p - y[i]) * wt[i];
    for (int k = 0; k < l2; ++k) {
      g_loc[k] += X[k * l1 + i] * resid;
    }
  }

  *qf_out = res_acc;
  std::memcpy(grad_out, g_loc, sizeof(double) * static_cast<size_t>(l2));
}

}  // namespace glmbayes::cl_ref

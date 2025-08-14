// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed



__kernel void f2_gamma(
    __global const double* X,      // design matrix,     l1 × l2, column‐major
    __global const double* B,      // grid matrix b,     l2 × m1, column‐major
    __global const double* mu,     // prior mean vector, length = l2
    __global const double* P,      // prior precision,   l2 × l2, row‐major
    __global const double* alpha,  // offset vector,     length = l1
    __global const double* y,      // responses,         length = l1
    __global const double* wt,     // weights = shape,   length = l1
    __global double*       qf,     // out: f2 + f3 values, length = m1
    __global double*       xb,     // out: scale = μ/shape, size = m1 × l1 (row‐major per grid)
    const int              l1,    // # observations
    const int              l2,    // # predictors
    const int              m1     // # grid points
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) PRIOR: tmp = P × (B_j – mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            double diff = B[ell * m1 + j] - mu[ell];
            acc += P[k*l2 + ell] * diff;
        }
        tmp[k] = acc;
    }

    // 2) PRIOR QUAD FORM: 0.5 * (B_j – mu)'·tmp
    double q_acc = 0.0;
    for (int k = 0; k < l2; ++k) {
        double diff = B[k * m1 + j] - mu[k];
        q_acc += diff * tmp[k];
    }
    q_acc *= 0.5;

    // 3) DATA TERM: for each observation
    //    η_i = α[i] + X[i,·]·B_j
    //    μ_i = exp(η_i)
    //    scale = μ_i / wt[i]
    //    ll = dgamma(y[i], shape = wt[i], scale, give_log=1)
    //    q_acc -= ll
    int off = j * l1;  // row‐major per grid: xb[off + i]
    for (int i = 0; i < l1; ++i) {
        // compute η
        double eta = alpha[i];
        for (int k = 0; k < l2; ++k) {
            eta += X[k*l1 + i] * B[k*m1 + j];
        }

        // link and scale
        double mui   = exp(eta);
        double scale = mui / wt[i];
        xb[off + i]  = scale;

        // accumulate negative‐log‐density
        double ll = dgamma(y[i], wt[i], scale, /*give_log=*/1);
        q_acc    -= ll;
    }

    // 4) WRITE‐BACK
    qf[j] = q_acc;
}
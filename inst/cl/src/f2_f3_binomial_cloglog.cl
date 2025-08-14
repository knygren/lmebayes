// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed


#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define MAX_L2 128

// f2 + f3 for Binomial–cloglog, single‐pass: prior + data‐term + gradient
__kernel void f2_f3_binomial_cloglog(
    __global const double* X,      // design matrix,    l1 × l2, column‐major
    __global const double* B,      // grid points,      m1 × l2, row‐major per grid
    __global const double* mu,     // prior mean,       length = l2
    __global const double* P,      // prior precision,  l2 × l2, row‐major
    __global const double* alpha,  // offsets,          length = l1
    __global const double* y,      // successes,        length = l1
    __global const double* wt,     // trials/weights,   length = l1
    __global double*       qf,     // out: neg-log-posterior, length = m1
    __global double*       xb,     // out: cloglog(p),       size = m1 × l1
    __global double*       grad,   // out: ∂(neg-log-post)/∂B, size = m1 × l2 (col-major)
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) Prior term: tmp[k] = [P × (B_j – mu)]_k
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    // 2) Quadratic form: 0.5 * (B_j – mu)' P (B_j – mu)
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double d_k = B[j*l2 + k] - mu[k];
        qsum += d_k * tmp[k];
    }
    double res_acc = 0.5 * qsum;

    // 3) Gradient accumulator starts with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) Data term: loop over observations
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // linear predictor η_i = α[i] + X[i,·]·B_j
        double eta = alpha[i];
        for (int k = 0; k < l2; ++k) {
            eta += X[k*l1 + i] * B[j*l2 + k];
        }

        // cloglog link and density factor
        double exp_eta    = exp(eta);
        double exp_neg    = exp(-exp_eta);
        double p1         = 1.0 - exp_neg;           // cloglog inverse
        double density    = exp(eta - exp_eta);      // derivative factor

        xb[base + i] = p1;

        // use dbinom for log-likelihood
        double ll = dbinom(y[i], wt[i], p1, /*give_log=*/1);
        res_acc -= ll;

        // gradient residual: ∂ℓ/∂η times wt
        // dℓ/dη = (y * density/p1) - ((wt - y) * density/(1-p1))
        double resid = ((y[i] * density / p1)
                      - ((wt[i] - y[i]) * density / exp_neg))
                      * wt[i];

        // accumulate gradient
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }

    // 5) Write back results
    qf[j] = res_acc;
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}
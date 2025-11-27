

#pragma once

// Required headers
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
#include <tbb/mutex.h>
static tbb::mutex f2_mutex;
#endif
#include <string>
#include "rng_utils.h"  // for safe_runif()
#include <atomic>
#include <memory>


using namespace Rcpp;
using namespace RcppParallel;

//-----------------------------------------------------------------------------
// rnnorm_reg_worker: parallel sampler with envelope logic
//-----------------------------------------------------------------------------
struct rnnorm_reg_worker : public RcppParallel::Worker {
  // --- Inputs ---
  int n;
  
  RVector<double>       y_r;       // observed counts
  RMatrix<double>       x_r;       // design matrix
  RMatrix<double>       mu_r;      // mode matrix
  RMatrix<double>       P_r;       // precision matrix
  RVector<double>       alpha_r;   // predictor offset
  RVector<double>       wt_r;      // observation weights
  
  // Envelope components as thread-safe handles (no copies)
  RVector<double> PLSD_r;
  RVector<double> LLconst_r;
  RMatrix<double> loglt_r;
  RMatrix<double> logrt_r;
  RMatrix<double> cbars_r;
  
  
  //   arma::vec             PLSD;      // slice density
  // arma::vec             LLconst;   // envelope constants
  // arma::mat             loglt;     // envelope lower bounds
  // arma::mat             logrt;     // envelope upper bounds
  // arma::mat             cbars;     // envelope centers
  
  CharacterVector       family;    // GLM family
  CharacterVector       link;      // link function
  int                   progbar;   // progress bar toggle
  
  // --- Outputs ---
  RMatrix<double>       out;       // accepted draws
  RVector<double>       draws;     // trial counts
  int                   ncol;      // dimensionality

  // --- Optional test controls ---
  // shared atomic flag: set to 1 by any thread if it hits the cap
  std::shared_ptr<std::atomic<int>> any_maxdraw_flag; // default nullptr (no reporting)
  int                   max_draws;                   // -1 => no per-index cap
  
  // --- Constructor ---
  rnnorm_reg_worker(
    int n_,
    const RVector<double>& y_r_,
    const RMatrix<double>& x_r_,
    const RMatrix<double>& mu_r_,
    const RMatrix<double>& P_r_,
    const RVector<double>& alpha_r_,
    const RVector<double>& wt_r_,
    
    const RcppParallel::RVector<double>& PLSD_r_,
    const RcppParallel::RVector<double>& LLconst_r_,
    const RcppParallel::RMatrix<double>& loglt_r_,
    const RcppParallel::RMatrix<double>& logrt_r_,
    const RcppParallel::RMatrix<double>& cbars_r_,
    
    // const arma::vec& PLSD_,
    // const arma::vec& LLconst_,
    // const arma::mat& loglt_,
    // const arma::mat& logrt_,
    // const arma::mat& cbars_,
    const CharacterVector& family_,
    const CharacterVector& link_,
    int progbar_,
    RMatrix<double>& out_,
    RVector<double>& draws_,
    std::shared_ptr<std::atomic<int>> any_maxdraw_flag_ = nullptr, // optional shared flag
    int max_draws_ = -1                                              // optional per-index cap
  )
    : n(n_),
      y_r(y_r_), x_r(x_r_), mu_r(mu_r_), P_r(P_r_),
      alpha_r(alpha_r_), wt_r(wt_r_),
      PLSD_r(PLSD_r_), LLconst_r(LLconst_r_),
      loglt_r(loglt_r_), logrt_r(logrt_r_), cbars_r(cbars_r_),
      // PLSD(PLSD_), LLconst(LLconst_),
      // loglt(loglt_), logrt(logrt_), cbars(cbars_),
      family(family_), link(link_), progbar(progbar_),
      out(out_), draws(draws_), ncol(out_.ncol())
    , any_maxdraw_flag(any_maxdraw_flag_),
      max_draws(max_draws_)
    
      
  {}
  
  // --- Parallel Loop ---
  void operator()(std::size_t begin, std::size_t end);
};



//-----------------------------------------------------------------------------
// rindep_norm_gamma_worker: parallel Normal–Gamma simulation with envelope
//-----------------------------------------------------------------------------
struct rindep_norm_gamma_worker : public RcppParallel::Worker {
  // --- Inputs ---
  int n;
  
  // Likelihood inputs (thread-safe views)
  RcppParallel::RVector<double>       y_r;
  RcppParallel::RMatrix<double>       x_r;
  RcppParallel::RMatrix<double>       mu_r;
  RcppParallel::RMatrix<double>       P_r;
  RcppParallel::RVector<double>       alpha_r;
  RcppParallel::RVector<double>       wt_r;
  
  // Envelope components
  RcppParallel::RMatrix<double>       cbars_r;
  RcppParallel::RVector<double>       PLSD_r;
  RcppParallel::RMatrix<double>       loglt_r;
  RcppParallel::RMatrix<double>       logrt_r;
  
  // UB vectors
  RcppParallel::RVector<double>       lg_prob_factor_r;
  RcppParallel::RVector<double>       UB2min_r;
  
  // Scalars
  double shape3, rate2, disp_upper, disp_lower, RSS_Min;
  double max_New_LL_UB, max_LL_log_disp, lm_log1, lm_log2, lmc1, lmc2;
  
  // Cache (precomputed upstream)
  RcppParallel::RMatrix<double>       Pmat_r;
  RcppParallel::RMatrix<double>       Pmu_r;
  RcppParallel::RVector<double>       base_B0_r;
  RcppParallel::RMatrix<double>       base_A_r;
  
  // --- Outputs ---
  RcppParallel::RMatrix<double>       beta_out_r;   // n × l1
  RcppParallel::RVector<double>       disp_out_r;   // length n
  RcppParallel::RVector<double>       iters_out_r;  // length n
  RcppParallel::RVector<double>       weight_out_r; // length n
  
  // --- Constructor ---
  rindep_norm_gamma_worker(
    int n_,
    const RcppParallel::RVector<double>& y_r_,
    const RcppParallel::RMatrix<double>& x_r_,
    const RcppParallel::RMatrix<double>& mu_r_,
    const RcppParallel::RMatrix<double>& P_r_,
    const RcppParallel::RVector<double>& alpha_r_,
    const RcppParallel::RVector<double>& wt_r_,
    const RcppParallel::RMatrix<double>& cbars_r_,
    const RcppParallel::RVector<double>& PLSD_r_,
    const RcppParallel::RMatrix<double>& loglt_r_,
    const RcppParallel::RMatrix<double>& logrt_r_,
    const RcppParallel::RVector<double>& lg_prob_factor_r_,
    const RcppParallel::RVector<double>& UB2min_r_,
    double shape3_, double rate2_,
    double disp_upper_, double disp_lower_,
    double RSS_Min_,
    double max_New_LL_UB_, double max_LL_log_disp_,
    double lm_log1_, double lm_log2_,
    double lmc1_, double lmc2_,
    const RcppParallel::RMatrix<double>& Pmat_r_,
    const RcppParallel::RMatrix<double>& Pmu_r_,
    const RcppParallel::RVector<double>& base_B0_r_,
    const RcppParallel::RMatrix<double>& base_A_r_,
    RcppParallel::RMatrix<double>& beta_out_r_,
    RcppParallel::RVector<double>& disp_out_r_,
    RcppParallel::RVector<double>& iters_out_r_,
    RcppParallel::RVector<double>& weight_out_r_)
    : n(n_),
      y_r(y_r_), x_r(x_r_), mu_r(mu_r_), P_r(P_r_), alpha_r(alpha_r_), wt_r(wt_r_),
      cbars_r(cbars_r_), PLSD_r(PLSD_r_), loglt_r(loglt_r_), logrt_r(logrt_r_),
      lg_prob_factor_r(lg_prob_factor_r_), UB2min_r(UB2min_r_),
      shape3(shape3_), rate2(rate2_), disp_upper(disp_upper_), disp_lower(disp_lower_),
      RSS_Min(RSS_Min_), max_New_LL_UB(max_New_LL_UB_), max_LL_log_disp(max_LL_log_disp_),
      lm_log1(lm_log1_), lm_log2(lm_log2_), lmc1(lmc1_), lmc2(lmc2_),
      Pmat_r(Pmat_r_), Pmu_r(Pmu_r_), base_B0_r(base_B0_r_), base_A_r(base_A_r_),
      beta_out_r(beta_out_r_), disp_out_r(disp_out_r_), iters_out_r(iters_out_r_), weight_out_r(weight_out_r_) {}
  
  // --- Parallel Loop ---
  void operator()(std::size_t begin, std::size_t end);
};


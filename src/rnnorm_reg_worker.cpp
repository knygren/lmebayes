#include "rnnorm_reg_worker.h"

// Only add more includes if strictly needed
#include <cmath>         // for std::log or std::exp if used
#include "famfuncs.h"
#include "Set_Grid.h"



//#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
//static std::mutex f2_mutex;
//#endif

// mutex to protect Rcpp calls


  // operator() implements the parallel loop
  void rnnorm_reg_worker::operator()(std::size_t begin, std::size_t end) {
    // Per-call lightweight views; no allocations of K×p
    arma::vec  PLSD   (PLSD_r.begin(),    PLSD_r.length(),               false, true);
    arma::vec  LLconst(LLconst_r.begin(), LLconst_r.length(),            false, true);
    arma::mat  loglt  (loglt_r.begin(),   loglt_r.nrow(), loglt_r.ncol(), false, true);
    arma::mat  logrt  (logrt_r.begin(),   logrt_r.nrow(), logrt_r.ncol(), false, true);
    arma::mat  cbars  (cbars_r.begin(),   cbars_r.nrow(), cbars_r.ncol(), false, true);
    
    
    
    // Create Armadillo views directly from RMatrix/RVector memory
    arma::vec y2(y_r.begin(), y_r.length(), false);
    arma::vec alpha2(alpha_r.begin(), alpha_r.length(), false);
    arma::vec wt2(wt_r.begin(), wt_r.length(), false);
    
    arma::mat x2(x_r.begin(), x_r.nrow(), x_r.ncol(), false);
    arma::mat mu2(mu_r.begin(), mu_r.nrow(), mu_r.ncol(), false);
    arma::mat P2(P_r.begin(), P_r.nrow(), P_r.ncol(), false);
    
    
    // Precompute dimensions and envelope pieces
    int l1 = mu_r.nrow();
    
    
    // Convert family/link once per thread
    std::string fam2 = as<std::string>(family);
    std::string lnk2 = as<std::string>(link);
    
    
    // Thread‐local buffers and views
    std::vector<double> outtemp_buf(l1), cbartemp_buf(l1);
    arma::rowvec        outtemp2(outtemp_buf.data(),   l1, false);
    arma::rowvec        cbartemp2(cbartemp_buf.data(), l1, false);
    
    
    std::vector<double> btemp_buf(l1);
    arma::mat btemp2(btemp_buf.data(), l1, 1, false);
    RcppParallel::RMatrix<double> btemp_r(btemp_buf.data(), l1, 1); // optional: only if still needed
    
    
    arma::mat testtemp2(1, 1);  // Allocated directly on the heap
    arma::vec testll2(1, arma::fill::none);  // Uninitialized vector of size m1
    
    /////////////////////////////////////////////////////////
    
    
    //    Rcpp::Rcout << "1.0 Launching Worker: " << begin << std::endl;
    
    // Main loop over indices
    for (std::size_t i = begin; i < end; ++i) {

      draws[i] = 1.0;  
      
      
      //      Rcpp::Rcout << "i=" << i  << "\n";
      
      double a1 = 0.0;
      
      while (a1 == 0.0) {
        
        
        // 1) slice selection
        //double U  = R::runif(0.0, 1.0)
        double U = safe_runif();
        double a2 = 0.0;
        int    J  = 0;
        while (a2 == 0.0) {
          if (U <= PLSD[J]) {
            //if (U <= PLSD2[J]) {
            a2 = 1.0;
          } else {
            U -= PLSD[J];
            ++J;
          }
        }
        
        // 2) draw truncated‐normal candidates
        for (int j = 0; j < l1; ++j) {
          out(i, j) = ctrnorm_cpp(logrt(J, j),loglt(J, j),-cbars(J, j), 1.0 );
          //  out(i, j) = ctrnorm_cpp(logrt2(J, j),loglt2(J, j),-cbars2(J, j), 1.0 );
        }
        
        // 3) prepare for test
        for (int j = 0; j < l1; ++j) {
          outtemp_buf[j]  = out(i, j);
          cbartemp_buf[j] = cbars(J, j);
          //cbartemp_buf[j] = cbars2(J, j);
          
          
        }
        testtemp2 = outtemp2 * trans(cbartemp2);
        //      double U2 = R::runif(0.0, 1.0);
        
        double U2 = safe_runif();
        
        
        btemp2   = trans(outtemp2);
        
        // declare test here so it’s in scope below
        //double test;
        
        
        
        // 4) compute log‐lik and print test under lock
        {
         
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
         tbb::mutex::scoped_lock lock(f2_mutex);
#endif  
         
//          std::lock_guard<std::mutex> guard(f2_mutex);
          
          
          
          // compute testll for all families/links
          if (fam2 == "binomial") {
            //            if (lnk2 == "logit")      testll = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
            if (lnk2 == "logit") 
            {
              testll2 = f2_binomial_logit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
              
              //Rcpp::Rcout << "rmat version: " << testll2  << "\n";
              
              //testll2 = f2_binomial_logit_arma(btemp,y,x,mu,P,alpha,wt,0);
              
              //              Rcpp::Rcout << "arma version: " << testll2  << "\n";
              
              //              testll2 = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
              //              Rcpp::Rcout << "original version: " << testll2  << "\n";
              
            }
            //if (lnk2 == "logit")      testll = f2_binomial_logit_arma(btemp,y,x,mu,P,alpha,wt,0);
            
            //else if (lnk2 == "probit") testll = f2_binomial_probit(btemp,y,x,mu,P,alpha,wt,0);
            //                    else if (lnk2 == "probit") testll = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
            else if (lnk2 == "probit") 
            {
              
              testll2 = f2_binomial_probit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
              //                      Rcpp::Rcout << "rmat version: " << testll2  << "\n";
              
              //                      testll2 = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
              
              //                      Rcpp::Rcout << "arma version: " << testll2  << "\n";
              
              //                      testll2 = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
            }
            //                    else                       testll = f2_binomial_cloglog(btemp,y,x,mu,P,alpha,wt,0);
            //                    else                       testll = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);
            else    
            {
              
              testll2 = f2_binomial_cloglog_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
              //                                            Rcpp::Rcout << "rmat version: " << testll2  << "\n";
              
              //                      testll2 = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);
              
              //                                            Rcpp::Rcout << "arma version: " << testll2  << "\n";
              
            }
          }
          else if (fam2 == "quasibinomial") {
            //            if (lnk2 == "logit")      testll = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
            if (lnk2 == "logit")
              
            {
              testll2 = f2_binomial_logit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
              //              testll2 = f2_binomial_logit_arma(btemp,y,x,mu,P,alpha,wt,0);
              //              testll2 = f2_binomial_logit(btemp,y,x,mu,P,alpha,wt,0);
              
            }
            //            else if (lnk2 == "probit") testll = f2_binomial_probit(btemp,y,x,mu,P,alpha,wt,0);
            //            else if (lnk2 == "probit") testll = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
            else if (lnk2 == "probit") 
              
            {
              //              Rcout << "Enter f2"  << std::endl;
              
              testll2 = f2_binomial_probit_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
              //                        Rcpp::Rcout << "rmat version: " << testll2  << "\n";
              
              //              testll2 = f2_binomial_probit_arma(btemp,y,x,mu,P,alpha,wt,0);
              
              //                          Rcpp::Rcout << "arma version: " << testll2  << "\n";
              
              //            Rcout << "Exit f2"  << std::endl;
            }
            
            //            else                       testll = f2_binomial_cloglog(btemp,y,x,mu,P,alpha,wt,0);
            //            else                       testll = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);
            else
              
            {
              testll2 = f2_binomial_cloglog_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
              
              //              testll2 = f2_binomial_cloglog_arma(btemp,y,x,mu,P,alpha,wt,0);
              
            }
          }
          else if (fam2 == "poisson"   || fam2 == "quasipoisson") {
            
            //            testll = f2_poisson(btemp,y,x,mu,P,alpha,wt,0);
            //            testll = f2_poisson_arma(btemp,y,x,mu,P,alpha,wt,0);
            
            testll2 = f2_poisson_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            
            //            Rcpp::Rcout << "rmat version v2: " << testll2  << "\n";
            
            
            //            testll2 = f2_poisson_rmat(btemp,y,x,mu,P,alpha,wt,0);
            
            //            Rcpp::Rcout << "rmat version: " << testll2  << "\n";
            
            //            testll2 = f2_poisson_arma(btemp,y,x,mu,P,alpha,wt,0);
            
            //            Rcpp::Rcout << "arma version: " << testll2  << "\n";
            
            
            //            testll[0]=testll2[0];  
          }
          else if (fam2 == "Gamma") {
            //            testll = f2_gamma(btemp,y,x,mu,P,alpha,wt,0);
            //            testll = f2_gamma_arma(btemp,y,x,mu,P,alpha,wt,0);
            testll2 = f2_gamma_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            //                        Rcpp::Rcout << "rmat version v2: " << testll2  << "\n";
            //            testll2 = f2_gamma_arma(btemp,y,x,mu,P,alpha,wt,0);
            //                        Rcpp::Rcout << "arma version: " << testll2  << "\n";
          }
          else { // gaussian
            //            testll = f2_gaussian(btemp,y,x,mu,P,alpha,wt);
            //            testll = f2_gaussian_arma(btemp,y,x,mu,P,alpha,wt);
            
            //  Note: This Envelope based sampling method for the Gaussian
            //        is not currently used. May implement future option
            //        to use as this is of theoretica interest
            //        and can be used to validate upper bounds
            testll2 = f2_gaussian_rmat(btemp_r,y_r,x_r,mu_r,P_r,alpha_r,wt_r,0);
            //  Rcpp::Rcout << "rmat version: " << testll2  << "\n";
            //            testll2 = f2_gaussian_arma(btemp,y,x,mu,P,alpha,wt);
            //            Rcpp::Rcout << "arma version: " << testll2  << "\n";
            
          }
          
          // calculate and print the acceptance statistic
          //          double test = LLconst[J]+ testtemp2(0,0) - std::log(U2)- testll[0];
          double test = LLconst[J]+ testtemp2(0,0) - std::log(U2)- testll2[0];
          
          // 5) Accept/reject logic
 
       
        
          if (test >= 0.0) {
            
            a1 = 1.0;            // accept
            
          } else {
            
            // keep existing behavior: increment trial count
            draws[i] = draws[i] + 1.0;
            
            // effective cap: use max_draws when provided, otherwise use legacy 1000 for diagnostic
//            int cap = (max_draws >= 0) ? max_draws : 1000;
            
            // print exactly once when we hit the cap (use your existing mutex for thread-safety)
//            if (static_cast<int>(draws[i]) == cap) {
              if (max_draws>0 && static_cast<int>(draws[i]) >= max_draws) {
              tbb::mutex::scoped_lock lock(f2_mutex);
              Rcpp::Rcout << "[WARN] index=" << i << " reached draws=" << draws[i]
                          << " (cap=" << max_draws << ") — forcing a1=1.0 to avoid infinite loop\n";
              
              Rcpp::Rcout << "[DEBUG] Acceptance test breakdown:\n";
              Rcpp::Rcout << "  LLconst[" << J << "] = " << LLconst[J] << "\n";
              Rcpp::Rcout << "  testtemp2(0,0) = " << testtemp2(0,0) << "\n";
              Rcpp::Rcout << "  log(U2) = " << std::log(U2) << "\n";
              Rcpp::Rcout << "  testll2[0] = " << testll2[0] << "\n";
              Rcpp::Rcout << "  test = " << test << "\n";            
            
            
            }
            
            
            
            // when cap reached or exceeded, set the atomic flag (if provided) and force exit
//            if (static_cast<int>(draws[i]) >= cap) {
              if (max_draws>0 && static_cast<int>(draws[i]) >= max_draws) {
                if (any_maxdraw_flag) {
                any_maxdraw_flag->store(1, std::memory_order_relaxed);
              }
              a1 = 1.0;   // force acceptance / break out of while loop
            }
            
          }
          
          
          
        }
        
        
        
        
        
      } // while(a1)
    }   // for(i)
    //  Rcpp::Rcout << "Exiting Worker: " << end << std::endl;
    
  }     // operator()



  // --- rindep_norm_gamma_worker implementation ---
  void rindep_norm_gamma_worker::operator()(std::size_t begin, std::size_t end) {
    const int l2 = x_r.nrow();
    const int l1 = x_r.ncol();
    
    for (std::size_t i = begin; i < end; ++i) {
      // Thread-local buffers and views (no shared state)
      std::vector<double> out_buf(static_cast<std::size_t>(l1), 0.0);
      RcppParallel::RMatrix<double> out_row(out_buf.data(), 1,  l1);  // 1×l1
      RcppParallel::RMatrix<double> out_col(out_buf.data(), l1, 1);   // l1×1
      
      std::vector<double> theta_buf(static_cast<std::size_t>(l1), 0.0);
      RcppParallel::RMatrix<double> theta_row(theta_buf.data(), 1,  l1); // 1×l1
      RcppParallel::RMatrix<double> theta_col(theta_buf.data(), l1, 1);  // l1×1
      
      std::vector<double> cbars_col_buf(static_cast<std::size_t>(l1), 0.0);
      RcppParallel::RMatrix<double> cbars_small_col(cbars_col_buf.data(), l1, 1); // l1×1
      
      
      // Scaled weights: classical logic requires wt2 = wt / dispersion before likelihood
   //   Rcpp::NumericVector wt2_nv(l2);                  // thread-local
    //  RcppParallel::RVector<double> wt2_r_old(wt2_nv);     // matches f2_gaussian_rmat signature
      
      
      std::vector<double> wt2_buf(static_cast<std::size_t>(l2), 0.0);
      // Wrap as a 1‑column matrix (l1 × 1)
      RcppParallel::RMatrix<double> wt2_r(wt2_buf.data(), l2, 1);

            
      iters_out_r[i]  = 1.0;
      weight_out_r[i] = 1.0;
      
      int accept = 0;
      

      
      while (accept == 0) {
        // 1) Slice/component selection via PLSD
        double U = safe_runif();
        int J_idx = 0;
        double U_left = U;
        while (true) {
          if (U_left <= PLSD_r[J_idx]) break;
          U_left -= PLSD_r[J_idx];
          ++J_idx;
        }
        
        // 2) Draw truncated-normal beta row
        for (int j = 0; j < l1; ++j) {
          out_row(0, j) = ctrnorm_cpp(
            logrt_r(J_idx, j),
            loglt_r(J_idx, j),
            -cbars_r(J_idx, j),
            1.0
          );
        }
        
        // 3) Draw dispersion
        double dispersion = r_invgamma_safe(shape3, rate2, disp_upper, disp_lower);
        
        // 4) Solve theta (strict row-only)
        for (int j = 0; j < l1; ++j) cbars_small_col(j, 0) = cbars_r(J_idx, j);
        
        // RcppParallel::RMatrix<double> theta_sol_r =
        //   Inv_f3_with_disp_rmat(Pmat_r, Pmu_r, base_B0_r, base_A_r,
        //                         dispersion, cbars_small_col);
        
        arma::mat theta_sol = Inv_f3_with_disp_rmat_v2(Pmat_r, Pmu_r, base_B0_r, base_A_r,
                                                       dispersion, cbars_small_col);
        
        
        //for (int j = 0; j < l1; ++j) theta_row(0, j) = theta_sol_r(0, j);
        
        for (int j = 0; j < l1; ++j) {
          theta_row(0, j) = theta_sol(0, j);
        }
        
        // 5) Scale weights before likelihood
        for (int r = 0; r < l2; ++r) {
          wt2_r(r, 0) = wt_r[r] / dispersion;
        }
        
        

        
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
        tbb::mutex::scoped_lock lock(f2_mutex);
#endif  
        
        // Rcout << "Entering f2_gaussian_rmat_mat 1"  << std::endl;
        
        // 6) Likelihood calls (column views, pre-scaled weights)
        double LL_New2_scalar =
          -f2_gaussian_rmat_mat(theta_col, y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];

          // Rcout << "Entering f2_gaussian_rmat_mat 2"  << std::endl;
          
          double LL_Test_scalar =
          -f2_gaussian_rmat_mat(out_col,   y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];
          

          
          // 7) Upper bounds
          double U2     = safe_runif();
          double log_U2 = std::log(U2);
          
          double UB1 = LL_New2_scalar;
          for (int j = 0; j < l1; ++j)
            UB1 -= cbars_r(J_idx, j) * (out_row(0, j) - theta_row(0, j));
          
          double quad_sum = 0.0;
          for (int r = 0; r < l2; ++r) {
            double x_theta = 0.0;
            for (int c = 0; c < l1; ++c) x_theta += x_r(r, c) * theta_row(0, c);
            double resid  = (y_r[r] - alpha_r[r] - x_theta);
            double scaled = resid * std::sqrt(wt_r[r]);
            quad_sum += scaled * scaled;
                      }
          double UB2 = 0.5 * (1.0 / dispersion) * (quad_sum - RSS_Min);
          UB2 -= UB2min_r[J_idx];
          
          double theta_P_theta = 0.0;
          for (int r = 0; r < l1; ++r) {
            double acc = 0.0;
            for (int c = 0; c < l1; ++c) acc += P_r(r, c) * theta_row(0, c);
            theta_P_theta += theta_row(0, r) * acc;
          }
          double c_theta = 0.0;
          for (int j = 0; j < l1; ++j) c_theta += cbars_r(J_idx, j) * theta_row(0, j);
          double New_LL_J = -0.5 * theta_P_theta + c_theta;
          
          double UB3A = lg_prob_factor_r[J_idx] + lmc1 + lmc2 * dispersion - New_LL_J;
          double New_LL_log_disp = lm_log1 + lm_log2 * std::log(dispersion);
          double UB3B = (max_New_LL_UB - max_LL_log_disp + New_LL_log_disp)
            - (lmc1 + lmc2 * dispersion);
          
          double test1 = (LL_Test_scalar - UB1);
          double test  = test1 - (UB2 + UB3A + UB3B);
          test = test - log_U2;
          
          
          // Rcout << "Entering Output assignment"  << std::endl;
          
          // 8) Record outputs and accept/reject
          disp_out_r[i] = dispersion;
          for (int j = 0; j < l1; ++j) beta_out_r(i, j) = out_row(0, j);
          
          if (test >= 0.0) {
            accept = 1;
          } else {
            iters_out_r[i] = iters_out_r[i] + 1;
          }
      } // while (accept == 0)
    }   // for i
  }
  

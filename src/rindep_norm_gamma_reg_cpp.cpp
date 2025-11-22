// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "Set_Grid.h"
#include <math.h>
#include "rng_utils.h"  // for safe_runif()

#include "nmath_local.h"
#include "dpq_local.h"


using namespace Rcpp;


void progress_bar3(double x, double N)
{
  // how wide you want the progress meter to be
  int totaldotz=40;
  double fraction = x / N;
  // part of the progressmeter that's already "full"
  int dotz = round(fraction * totaldotz);
  
  Rcpp::Rcout.precision(3);
  Rcout << "\r                                                                 " << std::flush ;
  Rcout << "\r" << std::flush ;
  Rcout << std::fixed << fraction*100 << std::flush ;
  Rcout << "% [" << std::flush ;
  int ii=0;
  for ( ; ii < dotz;ii++) {
    Rcout << "=" << std::flush ;
  }
  // remaining part (spaces)
  for ( ; ii < totaldotz;ii++) {
    Rcout << " " << std::flush ;
  }
  // and back to line begin 
  
  Rcout << "]" << std::flush ;
  
  // and back to line begin 
  
  Rcout << "\r" << std::flush ;
  
}


double p_inv_gamma(double dispersion,double shape,double rate){
  
  return(1- R::pgamma(1/dispersion,shape,1/rate,TRUE,FALSE));
}



double  q_inv_gamma(double p,double shape,double rate,double disp_upper,double disp_lower){
  double p_upp=p_inv_gamma(disp_upper,shape,rate);
  double p_low=p_inv_gamma(disp_lower,shape,rate);
  double p1=p_low+p*(p_upp-p_low);
  double p2=1-p1;
  return(1/ R::qgamma(p2,shape,1/rate,TRUE,FALSE));
}

double r_invgamma(double shape,double rate,double disp_upper,double disp_lower){
  double p= R::runif(0,1);
  return(q_inv_gamma(p,shape,rate,disp_upper,disp_lower));
}


// Safe inverse-gamma CDF using nmath/rmath pgamma
double p_inv_gamma_safe(double dispersion,
                        double shape,
                        double rate) {
  // For X ~ InvGamma(shape, rate), Y = 1/X ~ Gamma(shape, rate)
  // So P(X <= d) = P(Y >= 1/d) = 1 - F_Y(1/d)
  double y = 1.0 / dispersion;
  
  // Call the ported pgamma (not R::pgamma)
  // Arguments: x, shape, scale, lower_tail, log_p
  double Fy = pgamma_local(y, shape, 1.0 / rate, /*lower_tail=*/1, /*log_p=*/0);
  
  return 1.0 - Fy;
}


double q_inv_gamma_safe(double p,
                        double shape,
                        double rate,
                        double disp_upper,
                        double disp_lower) {
  // Compute probabilities at the bounds using safe pgamma
  double p_upp = p_inv_gamma_safe(disp_upper, shape, rate);
  double p_low = p_inv_gamma_safe(disp_lower, shape, rate);

  // Map uniform p into [p_low, p_upp]
  double p1 = p_low + p * (p_upp - p_low);
  double p2 = 1.0 - p1;

  // Invert via safe qgamma (ported from nmath/rmath)
  return 1.0 / qgamma_local(p2, shape, 1.0 / rate, /*lower_tail=*/1, /*log_p=*/0);
}



// 
// // Declaration (e.g. in a header if needed)
// // double r_invgamma_safe(double shape, double rate,
// //                        double disp_upper, double disp_lower);
// 
// // Definition (in your .cpp file)
double r_invgamma_safe(double shape,
                       double rate,
                       double disp_upper,
                       double disp_lower) {
  // draw uniform(0,1) from thread‑local RNG
  double p = safe_runif();

  // invert CDF at p to get inverse‑gamma draw
  // q_inv_gamma must be pure C++ math, no R calls
  return q_inv_gamma(p, shape, rate, disp_upper, disp_lower);
}



// [[Rcpp::export(".rindep_norm_gamma_reg_std_cpp")]]

Rcpp::List  rindep_norm_gamma_reg_std_cpp(int n,NumericVector y,NumericMatrix x,
                                             NumericMatrix mu, /// This is typically standardized to be a zero vector
                                             NumericMatrix P, /// Part of prior precision shifted to the likelihood
                                             NumericVector alpha,NumericVector wt,
                                             Function f2,Rcpp::List  Envelope,
                                             Rcpp::List  gamma_list,
                                             Rcpp::List  UB_list,
                                             Rcpp::CharacterVector   family,Rcpp::CharacterVector   link, bool progbar=true)
{
  
  // 1. Grab the base environment
  Rcpp::Environment base = Rcpp::Environment::base_env();
  
  // 2. Pull out the 'interactive' function
  Rcpp::Function interactive = base["interactive"];
  
  
  int l1 = mu.nrow();
  int l2 = x.nrow();
  
  
  // Get various inputs frm the provided lists
  
  double shape3 =gamma_list["shape3"];
  double rate2 =gamma_list["rate2"];
  double disp_upper =gamma_list["disp_upper"];
  double disp_lower =gamma_list["disp_lower"];
  double RSS_ML =UB_list["RSS_ML"];
  double max_New_LL_UB =UB_list["max_New_LL_UB"];
  double max_LL_log_disp =UB_list["max_LL_log_disp"];
  double lm_log1 =UB_list["lm_log1"];
  double lm_log2 =UB_list["lm_log2"];
  double lmc1 =UB_list["lmc1"];
  double lmc2 =UB_list["lmc2"];
  NumericVector lg_prob_factor =UB_list["lg_prob_factor"];
  NumericMatrix cbars=Envelope["cbars"];
  
  
  NumericVector iters_out(n);
  NumericVector disp_out(n);
  NumericVector weight_out(n);
  NumericMatrix beta_out(n,l1);
  double dispersion;
  NumericVector wt2(l1);
  
  
  arma::vec wt1b(wt.begin(), x.nrow());
  
  
  NumericMatrix cbarst(cbars.ncol(),cbars.nrow());
  NumericMatrix thetabars(cbars.nrow(),cbars.ncol());
  NumericMatrix thetabars_new(1,cbars.ncol());
  
  NumericVector New_LL(cbars.nrow());
  
  
  
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);
  
  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);
  cbarstb=trans(cbarsb);
  
  arma::vec y2(y.begin(),l2);
  arma::vec alpha2(alpha.begin(),l2);
  arma::mat x2(x.begin(),l2,l1);
  arma::mat P2(P.begin(),l1,l1);
  
  double UB1;
  double UB2;
  double UB3A;
  double UB3B;
  double New_LL_log_disp;
  
  int a1=0;
  double test1=0;
  double test=0;
  NumericVector J(n);
  NumericVector draws(n);
  NumericMatrix out(1,l1);
  double a2=0;
  double U=0;
  double U2=0;
  
  NumericVector PLSD=Envelope["PLSD"];
  NumericMatrix loglt=Envelope["loglt"];
  NumericMatrix logrt=Envelope["logrt"];
  
  double RSS_Min=UB_list["RSS_Min"];
  NumericVector UB2min=UB_list["UB2min"];
  
//  NumericVector ub2_min=;
  
  
  
  
  // Build cache once outside the loop
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  
  for(int i=0;i<n;i++){

    Rcpp::checkUserInterrupt();
    
//    if(progbar==1){
//      progress_bar3(i, n-1);
//      if(i==n-1) {Rcpp::Rcout << "" << std::endl;}
//    }
    
    // 3. Test progbar *and* interactive()



    
    a1=0;
    iters_out[i]=1;  
    while(a1==0){

          
      
      // Simulate from discrete distribution
      
      U=R::runif(0.0, 1.0);
      a2=0;
      J(0)=0;    
      while(a2==0){
        if(U<=PLSD(J(0))) a2=1;
        if(U>PLSD(J(0))){ 
          U=U-PLSD(J(0));
          J(0)=J(0)+1;
          
        }
      }
      

            
      // Simulate for beta
      
      for(int j=0;j<l1;j++){  out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);          }
      
      

      // Update this to make distribution contingent on component of the grid
      
      dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);
      
      
      
      wt2=wt/dispersion;
      NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
      
      // Compute Adjusted theta (accounting for changed dispersion) - New tangency points
    
      arma::mat theta2 = Inv_f3_with_disp(cache, dispersion, transpose(cbars_small));
      thetabarsb_new = theta2;
      

      // theta2 =Inv_f3_gaussian(transpose(cbars_small), y,x, mu, P, alpha, wt2);  
      // thetabarsb_new=theta2;
      

      // Recompoute LL at the new gradient point
      NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);  
      
    
      
      U2=R::runif(0.0, 1.0);
      
      double log_U2=log(U2);
      NumericVector J_out=J;
      NumericVector b_out=out(0,_);
      arma::rowvec b_out2(b_out.begin(),l1,false);
      NumericVector thetabars_temp=thetabars_new(0,_); // Changed
      
      arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
      NumericVector cbars_temp=cbars(J_out(0),_);
      arma::vec  cbars_temp2(cbars_temp.begin(), l1);
      
      
      
      NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);
      

      
      // Block 1: UB1 
      //   Same form as in fixed dispersion case but thetabar is a function of the dispersion
      //   So all components that include thetabar must now be bounded as well
      
      arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
      UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
      
      //Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
      
      
      arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
      
//      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
      UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_Min);
      
      // Subtract UB2min --> Should improve acceptance
      
      UB2=UB2-UB2min[J_out(0)];
      
      
      // Block 3: UB3A (adjusts because probabilities of components in grid are different from original grid)
      // Investigate whether changing probabilities of grid components for proposal
      // allows us to do away with this term and to thereby improve the acceptance rate
      
      // This is likely time consuming part
      

      
      for(int j=J_out(0);j<(J_out(0)+1);j++){
        thetabars_temp=thetabars_new(0,_); // Changed
        
        
        cbars_temp=cbars(j,_);
        arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
        arma::vec  cbars_temp2(cbars_temp.begin(), l1);
        
        New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
                                    +trans(cbars_temp2)*thetabars_temp2);
        
      }
      

      // Modified UB3A 
      
      UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
      
      // Block 4: UB3B  
      
      New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
      
      UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);
      

      
      test1=LL_Test[0]-UB1;
        
      test= test1-(UB2+UB3A+UB3B);  // Should be all negative 
      

      test = test - log_U2;
      
      disp_out[i] = dispersion;
      beta_out(i, _) = out(0, _);
      

      if(test>=0){
        

        
        a1=1;
        
      }
      else{
        iters_out[i]=iters_out[i]+1;
        }    
      

    }  
    
    
  }
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}




// Interleaved, line-by-line correspondence: each classical step commented, followed by the active rmat step.
// Ensures identical sequencing, orientation, and arithmetic — no extra intermediates.

// Fully revised, line-by-line mirrored implementation.
// Classical steps are commented immediately above the active rmat lines.
// Matrices/vectors, shapes, transposes, and copy semantics are identical.
// out is a 1×l1 matrix; theta_row is normalized to 1×l1 matrix before transpose;
// Likelihoods take l1×1 column inputs; acceptance test is a single line.

// Literal mirror of the classic and worker “column contract”
// - Draw row candidates
// - Transpose to strict l1×1 columns via shared buffers
// - Single wrap of columns passed to f2_gaussian_rmat
// - UB math identical to classic
// - Acceptance logic fixed (do not accept in the else branch)

// Thread-safe rindep_loop_rmat without any Rcpp::NumericVector/Matrix allocations.
// Key points:
// - Use std::vector<double> buffers + RcppParallel::RMatrix views (row + col) for beta and theta.
// - Enforce solver row-only output (1 × l1).
// - Pass wt_r directly to f2_gaussian_rmat; scale the returned LL by 1/dispersion to avoid building wt/dispersion.
// - UB math uses row views; likelihood uses column views.
// - No changes to f2_gaussian_rmat signature.

// Thread-safe rindep_loop_rmat
// - No Rcpp::NumericMatrix allocations for out/theta/cbars
// - Enforce solver row-only output (1 × l1)
// - Scale weights before likelihood call (wt2 = wt / dispersion)
// - Pass wt2 directly into f2_gaussian_rmat (no scaling of returned LL)

void rindep_loop_rmat(
    int n,
    // Inputs
    const RcppParallel::RVector<double>& y_r,
    const RcppParallel::RMatrix<double>& x_r,
    const RcppParallel::RMatrix<double>& mu_r,
    const RcppParallel::RMatrix<double>& P_r,
    const RcppParallel::RVector<double>& alpha_r,
    const RcppParallel::RVector<double>& wt_r,
    
    const RcppParallel::RMatrix<double>& cbars_r,
    const RcppParallel::RVector<double>& PLSD_r,
    const RcppParallel::RMatrix<double>& loglt_r,
    const RcppParallel::RMatrix<double>& logrt_r,
    
    const RcppParallel::RVector<double>& lg_prob_factor_r,
    const RcppParallel::RVector<double>& UB2min_r,
    
    double shape3,
    double rate2,
    double disp_upper,
    double disp_lower,
    double RSS_Min,
    double max_New_LL_UB,
    double max_LL_log_disp,
    double lm_log1,
    double lm_log2,
    double lmc1,
    double lmc2,
    
    // Cache (precomputed upstream)
    const RcppParallel::RMatrix<double>& Pmat_r,
    const RcppParallel::RMatrix<double>& Pmu_r,
    const RcppParallel::RVector<double>& base_B0_r,
    const RcppParallel::RMatrix<double>& base_A_r,
    
    // Outputs
    RcppParallel::RMatrix<double>& beta_out_r,   // n × l1
    RcppParallel::RVector<double>& disp_out_r,   // length n
    RcppParallel::RVector<double>& iters_out_r,  // length n
    RcppParallel::RVector<double>& weight_out_r  // length n
) {
  const int l2 = x_r.nrow();
  const int l1 = x_r.ncol();
  
  // Thread-local buffers and views (no Rcpp::NumericMatrix allocations).
  std::vector<double> out_buf(static_cast<std::size_t>(l1), 0.0);
  RcppParallel::RMatrix<double> out_row(out_buf.data(), 1,  l1);  // 1×l1
  RcppParallel::RMatrix<double> out_col(out_buf.data(), l1, 1);   // l1×1
  
  std::vector<double> theta_buf(static_cast<std::size_t>(l1), 0.0);
  RcppParallel::RMatrix<double> theta_row(theta_buf.data(), 1,  l1); // 1×l1
  RcppParallel::RMatrix<double> theta_col(theta_buf.data(), l1, 1);  // l1×1
  
  std::vector<double> cbars_col_buf(static_cast<std::size_t>(l1), 0.0);
  RcppParallel::RMatrix<double> cbars_small_col(cbars_col_buf.data(), l1, 1); // l1×1
  
  // Scaled weights: classical logic requires wt2 = wt / dispersion before likelihood
  Rcpp::NumericVector wt2_nv(l2);                  // thread-local
  RcppParallel::RVector<double> wt2_r(wt2_nv);     // matches f2_gaussian_rmat signature
  
  for (int i = 0; i < n; ++i) {
    
    Rcpp::checkUserInterrupt();  // allow user to break out
    
    iters_out_r[i] = 1;
    weight_out_r[i] = 1.0;
    
    int a1 = 0;
    while (a1 == 0) {
      // 1) Component selection via PLSD
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
      

      
      // 4) Solve theta (strict row-only contract)
      for (int j = 0; j < l1; ++j) cbars_small_col(j, 0) = cbars_r(J_idx, j);
      
      RcppParallel::RMatrix<double> theta_sol_r =
        Inv_f3_with_disp_rmat(Pmat_r, Pmu_r, base_B0_r, base_A_r,
                              dispersion, cbars_small_col);
      
      // Signedness-safe check (cast l1 to std::size_t for compare)
      if (!(theta_sol_r.nrow() == 1 &&
          static_cast<std::size_t>(theta_sol_r.ncol()) == static_cast<std::size_t>(l1))) {
        Rcpp::stop(std::string("Inv_f3_with_disp_rmat must return 1×l1 row; got ") +
          std::to_string(theta_sol_r.nrow()) + "x" +
          std::to_string(theta_sol_r.ncol()));
      }
      for (int j = 0; j < l1; ++j) theta_row(0, j) = theta_sol_r(0, j);
      
      // 5) Scale weights before likelihood (classical behavior)
      for (int r = 0; r < l2; ++r) wt2_r[r] = wt_r[r] / dispersion;
      
      // 6) Likelihood calls (column views, pre-scaled weights)
      
      

      double LL_New2_scalar =
        -f2_gaussian_rmat(theta_col, y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];
        double LL_Test_scalar =
        -f2_gaussian_rmat(out_col,   y_r, x_r, mu_r, P_r, alpha_r, wt2_r, 0)[0];
        
        // 7) Upper bounds (row math)
        double U2 = safe_runif();
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
  
  
  double test = test1 - (UB2 + UB3A + UB3B);
  

        test = test - log_U2;
        
      
        // 8) Record outputs and accept/reject
        disp_out_r[i] = dispersion;
        for (int j = 0; j < l1; ++j) beta_out_r(i, j) = out_row(0, j);
        
        if (test >= 0.0) {
          a1 = 1;
        } else {
          
          iters_out_r[i] = iters_out_r[i] + 1;
        }
    } // while (a1 == 0)
  } // for i
} // end rindep_loop_rmat  



    
// Classic loop implementation: consumes pre-extracted inputs.
// Calls f2_gaussian(...) directly (assumed defined elsewhere).
void rindep_loop_classic(
    int n,
    // Rcpp originals for f2_gaussian
    const Rcpp::NumericVector& y_nv,
    const Rcpp::NumericMatrix& x_nm,
    const Rcpp::NumericMatrix& mu_nm,
    const Rcpp::NumericMatrix& P_nm,
    const Rcpp::NumericVector& alpha_nv,
    const Rcpp::NumericVector& wt_nv,
    
    // Envelope matrices/vectors
    Rcpp::NumericMatrix& cbars,
    Rcpp::NumericVector& PLSD,
    Rcpp::NumericMatrix& loglt,
    Rcpp::NumericMatrix& logrt,
    
    // UB vectors
    Rcpp::NumericVector& lg_prob_factor,
    Rcpp::NumericVector& UB2min,
    
    // Scalar constants
    double shape3,
    double rate2,
    double disp_upper,
    double disp_lower,
    double RSS_Min,
    double max_New_LL_UB,
    double max_LL_log_disp,
    double lm_log1,
    double lm_log2,
    double lmc1,
    double lmc2,
    
    // Precomputed cache
    Rcpp::List& cache,
    
    // Armadillo views for UB math
    arma::vec& y2,          // length l2
    arma::vec& alpha2,      // length l2
    arma::mat& x2,          // l2 × l1
    arma::mat& P2,          // l1 × l1
    arma::vec& sqrt_wt1b,   // length l2: sqrt(wt)
    
    // Outputs
    Rcpp::NumericMatrix& beta_out,   // n × l1
    Rcpp::NumericVector& disp_out,   // length n
    Rcpp::NumericVector& iters_out,  // length n
    Rcpp::NumericVector& weight_out  // length n
) {

  const int l1 = x2.n_cols;
  const int l2 = x2.n_rows;
  
  
  // Declare storage before the for-loop
  int J_idx_first = -1;
  

  for (int i = 0; i < n; ++i) {
    
    
    int a1 = 0;
    iters_out[i] = 1;

    while (a1 == 0) {

      // Draw component index J via PLSD
      double U = safe_runif();
      int J_idx = 0;
      double U_left = U;
      while (true) {
        if (U_left <= PLSD[J_idx]) break;
        U_left -= PLSD[J_idx];
        ++J_idx;
      }

      // Simulate beta row
      Rcpp::NumericMatrix out(1, l1);
      for (int j = 0; j < l1; ++j) {
        out(0, j) = ctrnorm_cpp(logrt(J_idx, j), loglt(J_idx, j), -cbars(J_idx, j), 1.0);
      }

      // Dispersion draw
      double dispersion = r_invgamma_safe(shape3, rate2, disp_upper, disp_lower);


  
      
      // Compute theta row
      Rcpp::NumericMatrix cbars_small = cbars(Rcpp::Range(J_idx, J_idx),
                                              Rcpp::Range(0, cbars.ncol() - 1));
      arma::mat theta2 = Inv_f3_with_disp(cache, dispersion, Rcpp::transpose(cbars_small));


      
      Rcpp::NumericMatrix thetabars_new(1, l1);
      // Fill using theta2 exactly as returned (must be 1 × l1 row)
      if (!(theta2.n_rows == 1 && theta2.n_cols == l1)) {
        Rcpp::stop("Inv_f3_with_disp must return 1×l1 row; got " +
          std::to_string(theta2.n_rows) + "x" + std::to_string(theta2.n_cols) +
          ", expected 1×l1 (l1=" + std::to_string(l1) + ")");
      }
      
      // Deterministic copy into thetabars_new row
      for (int j = 0; j < l1; ++j) thetabars_new(0, j) = theta2(0, j);
      
      
      // Likelihoods (calls f2_gaussian directly with Rcpp inputs)
      Rcpp::NumericVector wt2(l2);
      for (int r = 0; r < l2; ++r) wt2[r] = wt_nv[r] / dispersion;


      Rcpp::NumericVector LL_New2 = -f2_gaussian(Rcpp::transpose(thetabars_new),
                                                 y_nv, x_nm, mu_nm, P_nm, alpha_nv, wt2);
      
      
      

      Rcpp::NumericVector LL_Test = -f2_gaussian(Rcpp::transpose(out),
                                                 y_nv, x_nm, mu_nm, P_nm, alpha_nv, wt2);

      double U2 = safe_runif();
      double log_U2 = std::log(U2);

      // UB1
      arma::rowvec b_out2(out.begin(), l1, false);
      arma::vec    theta_vec(thetabars_new.begin(), l1, false);
      Rcpp::NumericVector cbars_row = cbars(J_idx, Rcpp::_);
      arma::vec cbars_vec(cbars_row.begin(), l1, false);
      arma::colvec betadiff = b_out2.t() - theta_vec;
      double UB1 = LL_New2[0] - arma::as_scalar(cbars_vec.t() * betadiff);
      
      // UB2
      double quad_sum = 0.0;
      for (int r = 0; r < l2; ++r) {
        double x_theta = 0.0;
        for (int c = 0; c < l1; ++c) x_theta += x2(r, c) * theta_vec[c];
        double resid  = (y2[r] - alpha2[r] - x_theta);
        double scaled = resid * sqrt_wt1b[r];
        quad_sum += scaled * scaled;
      }
      double UB2 = 0.5 * (1.0 / dispersion) * (quad_sum - RSS_Min);
      UB2 -= UB2min[J_idx];
      
      // UB3A
      double theta_P_theta = arma::as_scalar(theta_vec.t() * P2 * theta_vec);
      double c_theta       = arma::as_scalar(cbars_vec.t() * theta_vec);
      double New_LL_J      = -0.5 * theta_P_theta + c_theta;
      double UB3A          = lg_prob_factor[J_idx] + lmc1 + lmc2 * dispersion - New_LL_J;
      
      // UB3B
      double New_LL_log_disp = lm_log1 + lm_log2 * std::log(dispersion);
      double UB3B = (max_New_LL_UB - max_LL_log_disp + New_LL_log_disp)
        - (lmc1 + lmc2 * dispersion);
      
      // Acceptance test
    //  double test1 = LL_Test[0] - UB1;
      
      double test1=LL_Test[0] - UB1;
      
      double test= test1-(UB2+UB3A+UB3B);  // Should be all negative 

      test = test - log_U2;
      
      
      // Record outputs
      disp_out[i] = dispersion;
      beta_out(i, Rcpp::_) = out(0, Rcpp::_);
      
      if (test >= 0.0) {
        a1 = 1;
      } else {
        iters_out[i] = iters_out[i] + 1;
      }
    } // end while
  }   // end for
  
  


}




// [[Rcpp::export(".rindep_norm_gamma_reg_std_parallel_cpp")]]

Rcpp::List rindep_norm_gamma_reg_std_parallel_cpp(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericMatrix mu,  // typically standardized to be a zero vector
    Rcpp::NumericMatrix P,   // part of prior precision shifted to the likelihood
    Rcpp::NumericVector alpha,
    Rcpp::NumericVector wt,
    Rcpp::Function f2,
    Rcpp::List Envelope,
    Rcpp::List gamma_list,
    Rcpp::List UB_list,
    Rcpp::CharacterVector family,
    Rcpp::CharacterVector link,
    bool progbar = true
) {
  // Base env (kept as-is)
  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function interactive = base["interactive"];
  
  const int l1 = mu.nrow();
  const int l2 = x.nrow();
  
  // Scalars from lists
  double shape3          = gamma_list["shape3"];
  double rate2           = gamma_list["rate2"];
  double disp_upper      = gamma_list["disp_upper"];
  double disp_lower      = gamma_list["disp_lower"];
  double RSS_ML          = UB_list["RSS_ML"];
  double max_New_LL_UB   = UB_list["max_New_LL_UB"];
  double max_LL_log_disp = UB_list["max_LL_log_disp"];
  double lm_log1         = UB_list["lm_log1"];
  double lm_log2         = UB_list["lm_log2"];
  double lmc1            = UB_list["lmc1"];
  double lmc2            = UB_list["lmc2"];
  
  Rcpp::NumericVector lg_prob_factor = UB_list["lg_prob_factor"];
  Rcpp::NumericMatrix cbars          = Envelope["cbars"];
  Rcpp::NumericVector PLSD           = Envelope["PLSD"];
  Rcpp::NumericMatrix loglt          = Envelope["loglt"];
  Rcpp::NumericMatrix logrt          = Envelope["logrt"];
  double RSS_Min                     = UB_list["RSS_Min"];
  Rcpp::NumericVector UB2min         = UB_list["UB2min"];
  
  // Outputs
  Rcpp::NumericVector iters_out(n);
  Rcpp::NumericVector disp_out(n);
  Rcpp::NumericVector weight_out(n);
  Rcpp::NumericMatrix beta_out(n, l1);
  
  // Locals
  double dispersion = 0.0;
  Rcpp::NumericVector wt2(l2); // length l2, if used elsewhere
  
  // Armadillo views as in original
  arma::vec wt1b(wt.begin(), x.nrow());
  Rcpp::NumericMatrix cbarst(cbars.ncol(), cbars.nrow());
  Rcpp::NumericMatrix thetabars(cbars.nrow(), cbars.ncol());
  Rcpp::NumericMatrix thetabars_new(1, cbars.ncol());
  Rcpp::NumericVector New_LL(cbars.nrow());
  
  arma::mat cbarsb(cbars.begin(), cbars.nrow(), cbars.ncol(), false);
  arma::mat cbarstb(cbarst.begin(), cbarst.nrow(), cbarst.ncol(), false);
  
  arma::mat thetabarsb(thetabars.begin(), thetabars.nrow(), thetabars.ncol(), false);
  arma::mat thetabarsb_new(thetabars_new.begin(), thetabars_new.nrow(), thetabars_new.ncol(), false);
  cbarstb = arma::trans(cbarsb);
  
  arma::vec y2(y.begin(), l2);
  arma::vec alpha2(alpha.begin(), l2);
  arma::mat x2(x.begin(), l2, l1);
  arma::mat P2(P.begin(), l1, l1);
  
  double UB1 = 0.0, UB2v = 0.0, UB3A = 0.0, UB3B = 0.0, New_LL_log_disp = 0.0;
  
  int a1 = 0;
  double test1 = 0.0;
  double test = 0.0;
  
  Rcpp::NumericVector J(n);
  Rcpp::NumericVector draws(n);
  Rcpp::NumericMatrix out(1, l1);
  double a2 = 0.0;
  double U  = 0.0;
  double U2 = 0.0;
  
  // Build cache once outside the loop
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  // Wrap outputs with RcppParallel views
  RcppParallel::RMatrix<double> beta_out_r(beta_out);
  RcppParallel::RVector<double> disp_out_r(disp_out);
  RcppParallel::RVector<double> iters_out_r(iters_out);
  RcppParallel::RVector<double> weight_out_r(weight_out);
  
  // Wrap inputs with RcppParallel views (used by rmat path)
  RcppParallel::RVector<double> y_r(y);
  RcppParallel::RMatrix<double> x_r(x);
  RcppParallel::RMatrix<double> mu_r(mu);
  RcppParallel::RMatrix<double> P_r(P);
  RcppParallel::RVector<double> alpha_r(alpha);
  RcppParallel::RVector<double> wt_r(wt);
  
  RcppParallel::RMatrix<double> cbars_r(cbars);
  RcppParallel::RVector<double> PLSD_r(PLSD);
  RcppParallel::RMatrix<double> loglt_r(loglt);
  RcppParallel::RMatrix<double> logrt_r(logrt);
  
  // Rcpp originals for classic loop
  const Rcpp::NumericVector& y_nv   = y;
  const Rcpp::NumericMatrix& x_nm   = x;
  const Rcpp::NumericMatrix& mu_nm  = mu;
  const Rcpp::NumericMatrix& P_nm   = P;
  const Rcpp::NumericVector& alpha_nv = alpha;
  const Rcpp::NumericVector& wt_nv    = wt;
  
  // sqrt(wt) for UB2 term
  arma::vec sqrt_wt1b(wt.begin(), l2, false);
  sqrt_wt1b = arma::sqrt(sqrt_wt1b);
  
  // Wrap cache components as RcppParallel views
  Rcpp::NumericMatrix Pmat_nm    = cache["Pmat"];
  Rcpp::NumericMatrix Pmu_nm     = cache["Pmu"];
  Rcpp::NumericVector base_B0_nv = cache["base_B0"];
  Rcpp::NumericMatrix base_A_nm  = cache["base_A"];
  
  RcppParallel::RMatrix<double> Pmat_r(Pmat_nm);
  RcppParallel::RMatrix<double> Pmu_r(Pmu_nm);
  RcppParallel::RVector<double> base_B0_r(base_B0_nv);
  RcppParallel::RMatrix<double> base_A_r(base_A_nm);
  
  // Wrap UB vectors for rmat call
  RcppParallel::RVector<double> lg_prob_factor_r(lg_prob_factor);
  RcppParallel::RVector<double> UB2min_r(UB2min);

  // Debug: entering rmat loop
//  Rcpp::Rcout << "[DEBUG] Entering rindep_loop_rmat\n";
  
    
  // Call the rmat loop to fill outputs first
  rindep_loop_rmat(
    n,
    // Likelihood inputs
    y_r, x_r, mu_r, P_r, alpha_r, wt_r,
    // Envelope components
    cbars_r, PLSD_r, loglt_r, logrt_r,
    // UB vectors
    lg_prob_factor_r,
    UB2min_r,
    // Scalars
    shape3, rate2, disp_upper, disp_lower,
    RSS_Min, max_New_LL_UB, max_LL_log_disp,
    lm_log1, lm_log2, lmc1, lmc2,
    // Cache
    Pmat_r, Pmu_r, base_B0_r, base_A_r,
    // Outputs
    beta_out_r, disp_out_r, iters_out_r, weight_out_r
  );
  
  // // --- Summaries after rmat loop
  // {
  //   // Column means for beta_out
  //   Rcpp::NumericVector beta_means(beta_out.ncol());
  //   for (int j = 0; j < beta_out.ncol(); ++j) {
  //     double sum = 0.0;
  //     for (int i = 0; i < beta_out.nrow(); ++i) sum += beta_out(i, j);
  //     beta_means[j] = sum / static_cast<double>(beta_out.nrow());
  //   }
  // 
  //   // mean(disp_out)
  //   double disp_sum = 0.0;
  //   for (int i = 0; i < disp_out.size(); ++i) disp_sum += disp_out[i];
  //   double disp_mean = disp_sum / static_cast<double>(disp_out.size());
  // 
  //   // mean(iters_out)
  //   double iters_sum = 0.0;
  //   for (int i = 0; i < iters_out.size(); ++i) iters_sum += iters_out[i];
  //   double iters_mean = iters_sum / static_cast<double>(iters_out.size());
  // 
  //   Rcpp::Rcout << "[SUMMARY][rmat] mean(beta_out): " << beta_means << "\n";
  //   Rcpp::Rcout << "[SUMMARY][rmat] mean(disp_out): " << disp_mean << "\n";
  //   Rcpp::Rcout << "[SUMMARY][rmat] mean(iters_out): " << iters_mean << "\n";
  // }
  // 
  // // --- Standard deviations after rmat loop
  // {
  //   Rcpp::NumericVector beta_sds(beta_out.ncol());
  //   for (int j = 0; j < beta_out.ncol(); ++j) {
  //     double mean_j = 0.0;
  //     for (int i = 0; i < beta_out.nrow(); ++i) mean_j += beta_out(i, j);
  //     mean_j /= static_cast<double>(beta_out.nrow());
  // 
  //     double var_j = 0.0;
  //     for (int i = 0; i < beta_out.nrow(); ++i) {
  //       double diff = beta_out(i, j) - mean_j;
  //       var_j += diff * diff;
  //     }
  //     beta_sds[j] = std::sqrt(var_j / static_cast<double>(beta_out.nrow() - 1));
  //   }
  //   Rcpp::Rcout << "[SUMMARY][rmat] sd(beta_out): " << beta_sds << "\n";
  // 
  //   double disp_mean = 0.0;
  //   for (int i = 0; i < disp_out.size(); ++i) disp_mean += disp_out[i];
  //   disp_mean /= static_cast<double>(disp_out.size());
  // 
  //   double disp_var = 0.0;
  //   for (int i = 0; i < disp_out.size(); ++i) {
  //     double diff = disp_out[i] - disp_mean;
  //     disp_var += diff * diff;
  //   }
  //   double disp_sd = std::sqrt(disp_var / static_cast<double>(disp_out.size() - 1));
  //   Rcpp::Rcout << "[SUMMARY][rmat] sd(disp_out): " << disp_sd << "\n";
  // }

  
  // // Debug: entering classic loop
  // Rcpp::Rcout << "[DEBUG] Entering rindep_loop_classic\n";
  // 
  // // Run the classic loop to fill outputs (for validation / comparison)
  // rindep_loop_classic(
  //   n,
  //   y_nv, x_nm, mu_nm, P_nm, alpha_nv, wt_nv,
  //   cbars, PLSD, loglt, logrt,
  //   lg_prob_factor, UB2min,
  //   shape3, rate2, disp_upper, disp_lower,
  //   RSS_Min, max_New_LL_UB, max_LL_log_disp,
  //   lm_log1, lm_log2, lmc1, lmc2,
  //   cache,
  //   y2, alpha2, x2, P2, sqrt_wt1b,
  //   beta_out, disp_out, iters_out, weight_out
  // );
  // 
  // // Summaries after classic loop
  // {
  //   Rcpp::NumericVector beta_means(beta_out.ncol());
  //   for (int j = 0; j < beta_out.ncol(); ++j) {
  //     double sum = 0.0;
  //     for (int i = 0; i < beta_out.nrow(); ++i) sum += beta_out(i, j);
  //     beta_means[j] = sum / static_cast<double>(beta_out.nrow());
  //   }
  //   
  //   double disp_sum = 0.0;
  //   for (int i = 0; i < disp_out.size(); ++i) disp_sum += disp_out[i];
  //   double disp_mean = disp_sum / static_cast<double>(disp_out.size());
  //   
  //   double iters_sum = 0.0;
  //   for (int i = 0; i < iters_out.size(); ++i) iters_sum += iters_out[i];
  //   double iters_mean = iters_sum / static_cast<double>(iters_out.size());
  //   
  //   Rcpp::Rcout << "[SUMMARY][classic] mean(beta_out): " << beta_means << "\n";
  //   Rcpp::Rcout << "[SUMMARY][classic] mean(disp_out): " << disp_mean << "\n";
  //   Rcpp::Rcout << "[SUMMARY][classic] mean(iters_out): " << iters_mean << "\n";
  // }
  // 
  // {
  //   Rcpp::NumericVector beta_sds(beta_out.ncol());
  //   for (int j = 0; j < beta_out.ncol(); ++j) {
  //     double mean_j = 0.0;
  //     for (int i = 0; i < beta_out.nrow(); ++i) mean_j += beta_out(i, j);
  //     mean_j /= static_cast<double>(beta_out.nrow());
  //     
  //     double var_j = 0.0;
  //     for (int i = 0; i < beta_out.nrow(); ++i) {
  //       double diff = beta_out(i, j) - mean_j;
  //       var_j += diff * diff;
  //     }
  //     beta_sds[j] = std::sqrt(var_j / static_cast<double>(beta_out.nrow() - 1));
  //   }
  //   Rcpp::Rcout << "[SUMMARY] sd(beta_out): " << beta_sds << "\n";
  //   
  //   double disp_mean = 0.0;
  //   for (int i = 0; i < disp_out.size(); ++i) disp_mean += disp_out[i];
  //   disp_mean /= static_cast<double>(disp_out.size());
  //   
  //   double disp_var = 0.0;
  //   for (int i = 0; i < disp_out.size(); ++i) {
  //     double diff = disp_out[i] - disp_mean;
  //     disp_var += diff * diff;
  //   }
  //   double disp_sd = std::sqrt(disp_var / static_cast<double>(disp_out.size() - 1));
  //   Rcpp::Rcout << "[SUMMARY] sd(disp_out): " << disp_sd << "\n";
  // }
  
  return Rcpp::List::create(
    Rcpp::Named("beta_out")   = beta_out,
    Rcpp::Named("disp_out")   = disp_out,
    Rcpp::Named("iters_out")  = iters_out,
    Rcpp::Named("weight_out") = weight_out
  );
}



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
      // Fill using theta2 exactly as returned (no reshaping/remapping)
      if (theta2.n_rows == l1 && theta2.n_cols == 1) {
        for (int j = 0; j < l1; ++j) thetabars_new(0, j) = theta2(j, 0);
      } else if (theta2.n_rows == 1 && theta2.n_cols == l1) {
        for (int j = 0; j < l1; ++j) thetabars_new(0, j) = theta2(0, j);
      } else {
        Rcpp::stop("theta2 dims mismatch: got " +
          std::to_string(theta2.n_rows) + "x" + std::to_string(theta2.n_cols) +
          ", expected 1xl1 or l1x1 (l1=" + std::to_string(l1) + ")");
      }
      
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
      double test1 = LL_Test[0] - UB1;
      double test  = test1 - (UB2 + UB3A + UB3B) - log_U2;
      
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

Rcpp::List  rindep_norm_gamma_reg_std_parallel_cpp(int n,NumericVector y,NumericMatrix x,
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
  
  
  
  // Allocate outputs
//  NumericMatrix beta_out(n, l1);
//  NumericVector disp_out(n);
//  NumericVector iters_out(n);
//  NumericVector weight_out(n);
  
  // Wrap with RcppParallel views
  RcppParallel::RMatrix<double> beta_out_r(beta_out);
  RcppParallel::RVector<double> disp_out_r(disp_out);
  RcppParallel::RVector<double> iters_out_r(iters_out);
  RcppParallel::RVector<double> weight_out_r(weight_out);
  
  // Wrap inputs you’ll need inside the loop
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

  //////////////////////  Loop Prep ///////////////////////////
  
  
  // --- Prepare Rcpp originals for f2_gaussian (used inside the classic loop)
  const Rcpp::NumericVector& y_nv  = y;
  const Rcpp::NumericMatrix& x_nm  = x;
  const Rcpp::NumericMatrix& mu_nm = mu;
  const Rcpp::NumericMatrix& P_nm  = P;
  const Rcpp::NumericVector& alpha_nv = alpha;
  const Rcpp::NumericVector& wt_nv    = wt;
  
  // --- Armadillo views for UB math (used inside the classic loop)
  // arma::vec y2(y.begin(), l2, false);
  // arma::vec alpha2(alpha.begin(), l2, false);
  // arma::mat x2(x.begin(), l2, l1, false);
  // arma::mat P2(P.begin(), l1, l1, false);
  
  // sqrt(wt) for UB2 term (keep as sqrt(wt); classic loop uses this)
  arma::vec sqrt_wt1b(wt.begin(), l2, false);
  sqrt_wt1b = arma::sqrt(sqrt_wt1b);
  
  // --- Envelope components already present above
  // Rcpp::NumericMatrix cbars = Envelope["cbars"];
  // Rcpp::NumericVector PLSD  = Envelope["PLSD"];
  // Rcpp::NumericMatrix loglt = Envelope["loglt"];
  // Rcpp::NumericMatrix logrt = Envelope["logrt"];
  
  // --- UB vectors already present above
  // double RSS_Min              = UB_list["RSS_Min"];
  // Rcpp::NumericVector UB2min  = UB_list["UB2min"];
  // Rcpp::NumericVector lg_prob_factor = UB_list["lg_prob_factor"];
  
  // --- Scalar constants already present above
  // double shape3, rate2, disp_upper, disp_lower,
  //        max_New_LL_UB, max_LL_log_disp, lm_log1, lm_log2, lmc1, lmc2;
  
  // --- Precomputed cache already built above
  // Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  
  ///////////////////// End of Loop Prep ///////////////////////
    
  ///////////////////// Call new Loop Function /////////////////
  
  // --- Before calling classic loop
  Rcpp::Rcout << "[DEBUG] Entering rindep_loop_classic\n";
  
  
  // --- Run the classic loop first to fill outputs
  rindep_loop_classic(
    n,
    // Rcpp originals for f2_gaussian
    y_nv, x_nm, mu_nm, P_nm, alpha_nv, wt_nv,
    // Envelope matrices/vectors
    cbars, PLSD, loglt, logrt,
    // UB vectors
    lg_prob_factor, UB2min,
    // Scalar constants
    shape3, rate2, disp_upper, disp_lower,
    RSS_Min, max_New_LL_UB, max_LL_log_disp,
    lm_log1, lm_log2, lmc1, lmc2,
    // Precomputed cache
    cache,
    // Armadillo views
    y2, alpha2, x2, P2, sqrt_wt1b,
    // Outputs (already allocated above)
    beta_out, disp_out, iters_out, weight_out
  );
  
  
//   Rcpp::Rcout << "[DEBUG] Exiting rindep_loop_classic\n";
//   
// 
//   // --- Summaries after classic loop
//   {
//     // beta_out column means
//     Rcpp::NumericVector beta_means(beta_out.ncol());
//     for (int j = 0; j < beta_out.ncol(); ++j) {
//       double sum = 0.0;
//       for (int i = 0; i < beta_out.nrow(); ++i) sum += beta_out(i, j);
//       beta_means[j] = sum / static_cast<double>(beta_out.nrow());
//     }
//     
//     // mean(disp_out)
//     double disp_sum = 0.0;
//     for (int i = 0; i < disp_out.size(); ++i) disp_sum += disp_out[i];
//     double disp_mean = disp_sum / static_cast<double>(disp_out.size());
//     
//     // mean(iters_out)
//     double iters_sum = 0.0;
//     for (int i = 0; i < iters_out.size(); ++i) iters_sum += iters_out[i];
//     double iters_mean = iters_sum / static_cast<double>(iters_out.size());
//     
//     Rcpp::Rcout << "[SUMMARY][classic] mean(beta_out): " << beta_means << "\n";
//     Rcpp::Rcout << "[SUMMARY][classic] mean(disp_out): " << disp_mean << "\n";
//     Rcpp::Rcout << "[SUMMARY][classic] mean(iters_out): " << iters_mean << "\n";
//   }
//   
//   {
//   
//   // --- Standard deviations for beta_out columns
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
//   // --- Standard deviation for disp_out
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
//   
  
//   // --- Before old loop
//   Rcpp::Rcout << "[DEBUG] Entering old loop\n";
//   
//   ///////////////////// End of Call to new Loop Function //////
//   
//   for(int i=0;i<n;i++){
//     
// 
//     
//     a1=0;
//     iters_out[i]=1;  
//     while(a1==0){
//       
//       // Simulate from discrete distribution
//       //U=R::runif(0.0, 1.0);
//       U = safe_runif();
//       
//       // [Needs replace] Use your thread‑safe RNG: safe_runif().
//       
//       a2=0;
//       J(0)=0;    
//       while(a2==0){
//         if(U<=PLSD(J(0))) a2=1;
//         if(U>PLSD(J(0))){ 
//           U=U-PLSD(J(0));
//           J(0)=J(0)+1;
//         }
//       }
//       
//       // Simulate for beta
//       for(int j=0;j<l1;j++){
//         out(0,j)=ctrnorm_cpp(logrt(J(0),j),loglt(J(0),j),-cbars(J(0),j),1.0);
//         // [OK] ctrnorm_cpp is already thread‑safe (uses safe_runif and safe_qnorm_logp).
//       }
//       
//       // Update this to make distribution contingent on component of the grid
//       //dispersion=r_invgamma(shape3,rate2,disp_upper,disp_lower);
//       dispersion=r_invgamma_safe(shape3,rate2,disp_upper,disp_lower);
//       // [Needs replace] Implement r_invgamma_safe (per‑thread engine) and call that.
//       
//       wt2=wt/dispersion;
//       // [OK] Pure elementwise arithmetic; make sure wt2 is a thread‑local buffer in parallel.
//       
//       NumericMatrix cbars_small = cbars( Range(J(0),J(0)) , Range(0,cbars.ncol()-1) );
//       // [Replace in step 2] Avoid NumericMatrix slicing inside the loop. Use RMatrix view or copy row into std::vector<double>.
//       
//       // Compute Adjusted theta (accounting for changed dispersion) - New tangency points
//       arma::mat theta2 = Inv_f3_with_disp(cache, dispersion, transpose(cbars_small));
//       // [Review/Replace] Provide Inv_f3_with_disp_rmat or a pure C++ version.
//       // Must be deterministic, no R calls, and avoid allocating/slicing per iteration.
//       
//       thetabarsb_new = theta2;
//       // [Replace in step 2] Avoid Armadillo assignment to shared R objects in a parallel context.
//       // Use thread‑local buffers and write to outputs via RMatrix/RVector.
//       
//       // Recompute LL at the new gradient point
//       NumericVector LL_New2=-f2_gaussian(transpose(thetabars_new),  y, x, mu, P, alpha, wt2);
//       
//       // Call the Armadillo-backed version without reversing the sign
//       arma::vec LL_New2_rmat = -f2_gaussian_rmat(
//         RcppParallel::RMatrix<double>(thetabars_new),
//         y_r, x_r, mu_r, P_r, alpha_r, wt_r, 0);
//       
//       // [Needs replace] Use f2_gaussian_rmat(...) (pure C++, RMatrix/RVector inputs).
//       // Also avoid transpose(thetabars_new); pass a std::vector<double> or pointer.
//       
// //      U2=R::runif(0.0, 1.0);
//       U2 = safe_runif();
//       
//       // [Needs replace] Use safe_runif().
//       
//       double log_U2=log(U2);
//       // [OK] std::log is thread‑safe.
//       
//       NumericVector J_out=J;
//       NumericVector b_out=out(0,_);
//       arma::rowvec b_out2(b_out.begin(),l1,false);
//       NumericVector thetabars_temp=thetabars_new(0,_); // Changed
//       arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
//       NumericVector cbars_temp=cbars(J_out(0),_);
//       arma::vec  cbars_temp2(cbars_temp.begin(), l1);
//       // [Replace in step 2] These create Armadillo/NumericVector views/slices.
//       // Use thread‑local std::vector<double> buffers populated via RMatrix/RVector reads.
//       
//       NumericVector LL_Test=-f2_gaussian(transpose(out),  y, x, mu, P, alpha, wt2);
//       // [Needs replace] Use f2_gaussian_rmat(beta_row, y_r, x_r, mu_r, P_r, alpha_r, wt_r, dispersion).
//       
//       arma::vec LL_Test_rmat = -f2_gaussian_rmat(
//         RcppParallel::RMatrix<double>(out),
//         y_r, x_r, mu_r, P_r, alpha_r, wt_r, 0
//       );
//       
//       
//       // Block 1: UB1 
//       arma::colvec betadiff=trans(b_out2)-thetabars_temp2;
//       UB1=LL_New2(0) -arma::as_scalar(trans(cbars_temp2)*betadiff);
//       // [Replace in step 2] Avoid Armadillo trans/as_scalar. Compute scalar sums via loops.
//       
//       // Block 2: UB2 [RSS Term bounded by shifting it to the gamma candidate]
//       arma::colvec yxbeta=(y2-alpha2-x2*thetabars_temp2)%sqrt(wt1b); 
//       // [Replace in step 2] Replace with explicit loops over rows: dot(x_row, theta) then residual*sqrt(wt).
//       
//       // UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_ML);
//       UB2=0.5*(1/dispersion)*(arma::as_scalar(trans(yxbeta)*yxbeta)-RSS_Min);
//       // [Replace in step 2] Compute the quadratic sum with a scalar accumulator; avoid Armadillo.
//       
//       // Subtract UB2min --> Should improve acceptance
//       UB2=UB2-UB2min[J_out(0)];
//       // [OK] Pure scalar/indexing.
//       
//       // Block 3: UB3A ...
//       for(int j=J_out(0);j<(J_out(0)+1);j++){
//         thetabars_temp=thetabars_new(0,_); // Changed
//         cbars_temp=cbars(j,_);
//         arma::vec  thetabars_temp2(thetabars_temp.begin(), l1);
//         arma::vec  cbars_temp2(cbars_temp.begin(), l1);
//         New_LL(j)=arma::as_scalar(-0.5*trans(thetabars_temp2)*P2*thetabars_temp2
//                                     +trans(cbars_temp2)*thetabars_temp2);
//         // [Replace in step 2] Compute theta^T P theta and c^T theta with loops (RMatrix P_r and std::vector<double> theta).
//       }
//       
//       // Modified UB3A 
//       UB3A= lg_prob_factor(J_out(0))+lmc1+lmc2*dispersion-New_LL(J_out(0));
//       // [OK] Pure scalar/indexing.
//       
//       // Block 4: UB3B  
//       New_LL_log_disp=lm_log1+lm_log2*log(dispersion);
//       // [OK] std::log.
//       
//       UB3B=(max_New_LL_UB-max_LL_log_disp+New_LL_log_disp)-(lmc1+lmc2*dispersion);
//       // [OK] Pure scalar arithmetic.
//       
//       test1=LL_Test[0]-UB1;
//       // [OK] Scalar arithmetic (post replacement of LL_Test computation).
//       
//       test= test1-(UB2+UB3A+UB3B);  // Should be all negative 
//       // [OK] Scalar arithmetic.
//       
//       test = test - log_U2;
//       // [OK] Scalar arithmetic.
//       
//       disp_out[i] = dispersion;
//       beta_out(i, _) = out(0, _);
//       // [Replace in step 2] In parallel, write via RVector/RMatrix views:
//       // disp_out_r[i] = dispersion; for (int j=0; j<l1; ++j) beta_out_r(i,j) = beta_row[j];
//       
//       if(test>=0){
//         a1=1;
//       }
//       else{
//         iters_out[i]=iters_out[i]+1;
//         // [OK] Scalar write; in parallel still fine if each thread writes its own i.
//       }
//     }  
//   }
//   
//   Rcpp::Rcout << "[DEBUG] Exiting old loop\n";
//   
//   
//   // --- Summaries after old loop
//   {
//     // beta_out column means
//     Rcpp::NumericVector beta_means(beta_out.ncol());
//     for (int j = 0; j < beta_out.ncol(); ++j) {
//       double sum = 0.0;
//       for (int i = 0; i < beta_out.nrow(); ++i) sum += beta_out(i, j);
//       beta_means[j] = sum / static_cast<double>(beta_out.nrow());
//     }
//     
//     // mean(disp_out)
//     double disp_sum = 0.0;
//     for (int i = 0; i < disp_out.size(); ++i) disp_sum += disp_out[i];
//     double disp_mean = disp_sum / static_cast<double>(disp_out.size());
//     
//     // mean(iters_out)
//     double iters_sum = 0.0;
//     for (int i = 0; i < iters_out.size(); ++i) iters_sum += iters_out[i];
//     double iters_mean = iters_sum / static_cast<double>(iters_out.size());
//     
//     Rcpp::Rcout << "[SUMMARY][old] mean(beta_out): " << beta_means << "\n";
//     Rcpp::Rcout << "[SUMMARY][old] mean(disp_out): " << disp_mean << "\n";
//     Rcpp::Rcout << "[SUMMARY][old] mean(iters_out): " << iters_mean << "\n";
//   }
//   
//   
//   {
//     
//     // --- Standard deviations for beta_out columns
//     Rcpp::NumericVector beta_sds(beta_out.ncol());
//     for (int j = 0; j < beta_out.ncol(); ++j) {
//       double mean_j = 0.0;
//       for (int i = 0; i < beta_out.nrow(); ++i) mean_j += beta_out(i, j);
//       mean_j /= static_cast<double>(beta_out.nrow());
//       
//       double var_j = 0.0;
//       for (int i = 0; i < beta_out.nrow(); ++i) {
//         double diff = beta_out(i, j) - mean_j;
//         var_j += diff * diff;
//       }
//       beta_sds[j] = std::sqrt(var_j / static_cast<double>(beta_out.nrow() - 1));
//     }
//     Rcpp::Rcout << "[SUMMARY] sd(beta_out): " << beta_sds << "\n";
//     
//     // --- Standard deviation for disp_out
//     double disp_mean = 0.0;
//     for (int i = 0; i < disp_out.size(); ++i) disp_mean += disp_out[i];
//     disp_mean /= static_cast<double>(disp_out.size());
//     
//     double disp_var = 0.0;
//     for (int i = 0; i < disp_out.size(); ++i) {
//       double diff = disp_out[i] - disp_mean;
//       disp_var += diff * diff;
//     }
//     double disp_sd = std::sqrt(disp_var / static_cast<double>(disp_out.size() - 1));
//     Rcpp::Rcout << "[SUMMARY] sd(disp_out): " << disp_sd << "\n";
//   } 
  
  // Temporarily just return non-sense constants equal to all 1
  
  return Rcpp::List::create(Rcpp::Named("beta_out")=beta_out,Rcpp::Named("disp_out")=disp_out,
                            Rcpp::Named("iters_out")=iters_out,Rcpp::Named("weight_out")=weight_out);  
  
  
  
}




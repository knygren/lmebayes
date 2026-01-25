// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "famfuncs.h"
#include "Envelopefuncs.h"
#include "kernel_wrappers.h"
#include <RcppParallel.h>
#include "openclPort.h"
#include "utils_timing.h"

using namespace Rcpp;
using namespace openclPort;
using namespace famfuncs;







// [[Rcpp::export(".EnvelopeBuild_Ind_Normal_Gamma")]]

List EnvelopeBuild_Ind_Normal_Gamma(NumericVector bStar,NumericMatrix A,
                                    NumericVector y, 
                                    NumericMatrix x,
                                    NumericMatrix mu,
                                    NumericMatrix P,
                                    NumericVector alpha,
                                    NumericVector wt,
                                    std::string family,
                                    std::string link,
                                    int Gridtype, 
                                    int n,
                                    int n_envopt,
                                    bool sortgrid,
                                    bool use_opencl    ,
                                    bool verbose       
){
  
  
  //  int progbar=0;
  
  int l1 = A.nrow(), k = A.ncol();
  arma::mat A2(A.begin(), l1, k, false);
  arma::colvec bStar_2(bStar.begin(), bStar.size(), false);
  
  
  NumericVector a_1(l1);
  arma::vec a_2(a_1.begin(), a_1.size(), false);
  
  NumericVector xx_1(3, 1.0);
  NumericVector xx_2=NumericVector::create(-1.0,0.0,1.0);
  NumericVector yy_1(2, 1.0);
  NumericVector yy_2=NumericVector::create(-0.5,0.5);
  NumericMatrix G1(3,l1);
  NumericMatrix Lint1(2,l1);
  arma::mat G1b(G1.begin(), 3, l1, false);
  arma::mat Lint(Lint1.begin(), 2, l1, false);
  
  arma::colvec xx_1b(xx_1.begin(), xx_1.size(), false);
  arma::colvec xx_2b(xx_2.begin(), xx_2.size(), false);
  arma::colvec yy_1b(yy_1.begin(), yy_1.size(), false);
  arma::colvec yy_2b(yy_2.begin(), yy_2.size(), false);
  List G2(a_1.size());
  List GIndex1(a_1.size());
  Rcpp::Function EnvelopeOpt("EnvelopeOpt");
  Rcpp::Function expGrid("expand.grid");
  Rcpp::Function asMat("as.matrix");
  Rcpp::Function EnvSort("EnvelopeSort");
  
  int i;  
  
  a_2=arma::diagvec(A2);
  arma::vec omega=(sqrt(2)-arma::exp(-1.20491-0.7321*sqrt(0.5+a_2)))/arma::sqrt(1+a_2);
  G1b=xx_1b*arma::trans(bStar_2)+xx_2b*arma::trans(omega);
  Lint=yy_1b*arma::trans(bStar_2)+yy_2b*arma::trans(omega);
  
  // Second row in G1b here is the posterior mode
  
  NumericVector gridindex(l1);
  
  if(Gridtype==2){
    gridindex=EnvelopeOpt(a_2,n);
  }
  
  NumericVector Temp1=G1( _, 0);
  double Temp2;
  
  // Should write a small note with logic behind types 1 and 2
  
  for(i=0;i<l1;i++){
    
    if(Gridtype==1){
      
      // For Gridtype==1, small 1+a[i]<=(2/sqrt(M_PI) yields grid over full line
      // Can check speed for simulation when Gridtype=1 vs. Gridtyp=2 or 3     
      
      if((1+a_2[i])<=(2/sqrt(M_PI))){ 
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
      }
      if((1+a_2[i])>(2/sqrt(M_PI))){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      }    
    }  
    if(Gridtype==2){
      if(gridindex[i]==1){
        Temp2=G1(1,i);
        G2[i]=NumericVector::create(Temp2);
        GIndex1[i]=NumericVector::create(4.0);
      }
      if(gridindex[i]==3){
        Temp1=G1(_,i);
        G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
        GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
      }
    }
    
    if(Gridtype==3){
      Temp1=G1(_,i);
      G2[i]=NumericVector::create(Temp1(0),Temp1(1),Temp1(2));
      GIndex1[i]=NumericVector::create(1.0,2.0,3.0);
    }
    
    if(Gridtype==4){
      Temp2=G1(1,i);
      G2[i]=NumericVector::create(Temp2);
      GIndex1[i]=NumericVector::create(4.0);
    }
    
    
    
  }
  
  NumericMatrix G3=asMat(expGrid(G2));
  NumericMatrix GIndex=asMat(expGrid(GIndex1));
  NumericMatrix G4(G3.ncol(),G3.nrow());
  int l2=GIndex.nrow();
  
  arma::mat G3b(G3.begin(), G3.nrow(), G3.ncol(), false);
  arma::mat G4b(G4.begin(), G4.nrow(), G4.ncol(), false);
  
  G4b=trans(G3b);
  
  NumericMatrix cbars(l2,l1);
  NumericMatrix cbars_slope(l2,l1);
  NumericMatrix Up(l2,l1);
  NumericMatrix Down(l2,l1);
  NumericMatrix logP(l2,2);
  NumericMatrix logU(l2,l1);
  NumericMatrix loglt(l2,l1);
  NumericMatrix logrt(l2,l1);
  NumericMatrix logct(l2,l1);
  
  NumericMatrix LLconst(l2,1);
  NumericVector NegLL(l2);    
  NumericVector NegLL_slope(l2);    
  NumericVector RSS_Out(l2);
  arma::mat cbars2(cbars.begin(), l2, l1, false); 
  arma::mat cbars3(cbars.begin(), l2, l1, false); 
  
  arma::mat cbars_slope2(cbars_slope.begin(), l2, l1, false); 
  arma::mat cbars_slope3(cbars_slope.begin(), l2, l1, false); 
  
  
  // Note: NegLL_2 only added to allow for QC printing of results 
  
  arma::colvec NegLL_2(NegLL.begin(), NegLL.size(), false);
  
  //    G4b.print("tangent points");
  
  //  Rcpp::Rcout << "Gridtype is :"  << Gridtype << std::endl;
  //  Rcpp::Rcout << "Number of Variables in model are :"  << l1 << std::endl;
  //  Rcpp::Rcout << "Number of points in Grid are :"  << l2 << std::endl;
  
  
  
  
  if(family=="gaussian" ){
    //Rcpp::Rcout << "Finding Values of Log-posteriors:" << std::endl;
    
    // Adjust the slope calculations to split into several terms:
    // (i) Terms from shifted "prior" that does not depend on the dispersion
    // (ii) Constant terms from the actual LL that do not depend on dispersion or beta
    // (iii) Term from the LL that depends on the dispersion but not beta
    // (iv) Term from the LL that depends on beta and the dispersion (scaled RSS)
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Starting EnvelopeEval (NegLL, cbars) at " << now_hms() << " <<<\n";
    }
    Timer t_eval1; if (verbose) t_eval1.begin();
    
    Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt, family, link, use_opencl, verbose);
    NegLL = eval_info["NegLL"];
    cbars2 = Rcpp::as<arma::mat>(eval_info["cbars"]);
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting EnvelopeEval (NegLL, cbars) at " << now_hms() << " <<<\n";
      print_completed("[EnvelopeBuild] EnvelopeEval (NegLL, cbars)", t_eval1);
    }
    
    
//    Rcpp::List eval_info = EnvelopeEval(G4, y, x, mu, P, alpha, wt,
//                                        family, link, use_opencl, verbose);
    
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Starting EnvelopeEval (slope variants) at " << now_hms() << " <<<\n";
    }
    Timer t_eval2; if (verbose) t_eval2.begin();
    
    Rcpp::List eval_info2 = EnvelopeEval(G4, y, x, mu, 0*P, alpha, wt, family, link, use_opencl, verbose);
    NegLL_slope  = eval_info2["NegLL"];
    cbars_slope2 = Rcpp::as<arma::mat>(eval_info2["cbars"]);
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting EnvelopeEval (slope variants) at " << now_hms() << " <<<\n";
      print_completed("[EnvelopeBuild] EnvelopeEval (slope variants)", t_eval2);
      Rcpp::Rcout << "[EnvelopeBuild] Finished assigning NegLL_slope and cbars_slope2\n";
    }
    
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Starting RSS evaluation at " << now_hms() << " <<<\n";
    }
    Timer t_rss; if (verbose) t_rss.begin();
    
    RSS_Out = RSS(y, x, G4, alpha, wt); // includes dispersion in weight
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting RSS evaluation at " << now_hms() << " <<<\n";
      print_completed("[EnvelopeBuild] RSS evaluation", t_rss);
    } 
  }
  
  
  //  Rcpp::Rcout << "Finished Log-posterior evaluations:" << std::endl;
  
  // Do a temporary correction here cbars3 should point to correct memory
  // See if this sets cbars
  
  cbars3=cbars2;
  cbars_slope3=cbars_slope2;

  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Entering Set_Grid_C2 at " << now_hms() << " <<<\n";
  }
  Timer t_setgrid; if (verbose) t_setgrid.begin();
  
  Set_Grid_C2(GIndex, cbars, Lint1, Down, Up, loglt, logrt, logct, logU, logP);
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting Set_Grid_C2 at " << now_hms() << " <<<\n";
    print_completed("[EnvelopeBuild] Set_Grid_C2", t_setgrid);
  }
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Entering Set_logP_C2 at " << now_hms() << " <<<\n";
  }
  Timer t_setlogp; if (verbose) t_setlogp.begin();
  
  setlogP_C2(logP, NegLL, cbars, G3, LLconst);
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting Set_logP_C2 at " << now_hms() << " <<<\n";
    print_completed("[EnvelopeBuild] Set_logP_C2", t_setlogp);
  }  
  
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Starting PLSD computation at " << now_hms() << " <<<\n";
  }
  Timer t_plsd; if (verbose) t_plsd.begin();
  
  NumericMatrix::Column logP2 = logP(_, 1);
  double maxlogP = max(logP2);
  NumericVector PLSD = exp(logP2 - maxlogP);
  double sumP = sum(PLSD);
  PLSD = PLSD / sumP;
  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeBuild] >>> Exiting PLSD computation at " << now_hms() << " <<<\n";
    print_completed("[EnvelopeBuild] PLSD computation", t_plsd);
  }
  
  
  
  // Add sorting step back later after modifying EnvSort function
  // Should accomodate ready List
  
  //  if(sortgrid==true){
  //    Rcpp::List outlist=EnvSort(l1,l2,GIndex,G3,cbars,logU,logrt,loglt,logP,LLconst,PLSD,a_1);
  //    return(outlist);
  //  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("GridIndex")=GIndex,
                            Rcpp::Named("thetabars")=G3,
                            Rcpp::Named("cbars")=cbars,
                            Rcpp::Named("cbars_slope")=cbars_slope,
                            Rcpp::Named("NegLL")=NegLL,
                            Rcpp::Named("NegLL_slope")=NegLL_slope,
                            Rcpp::Named("Lint1")=Lint1,
                            Rcpp::Named("RSS_Out")=RSS_Out,
                            Rcpp::Named("logU")=logU,
                            Rcpp::Named("logrt")=logrt,
                            Rcpp::Named("loglt")=loglt,
                            Rcpp::Named("LLconst")=LLconst,
                            Rcpp::Named("logP")=logP(_,0),
                            Rcpp::Named("PLSD")=PLSD,
                            Rcpp::Named("a1")=a_1
  );
  
  
}



NumericVector RSS(NumericVector y, NumericMatrix x,NumericMatrix b,NumericVector alpha,NumericVector wt)
{
  // Step 1: Set up dimensions
  
  int l1 = x.nrow(), l2 = x.ncol(); // Dimensions of x matrix (dims for y,alpha, and wt needs to be consistent) 
  int m1 = b.ncol();                // Number of columns for which output is needed
  
  // Step 2: Initialize b2temp and other Rcpp and arma objects used in calculations
  
  Rcpp::NumericMatrix b2temp(l2,1);
  Rcpp::NumericMatrix restemp(1,1);
  arma::mat y2(y.begin(), l1, 1, false);
  arma::mat x2(x.begin(), l1, l2, false); 
  arma::mat alpha2(alpha.begin(), l1, 1, false); 
  
  Rcpp::NumericVector xb(l1);
  arma::colvec xb2(xb.begin(),l1,false); // Reuse memory - update both below
  
  NumericVector sqrt_wt=sqrt(wt);
  arma::mat sqrt_wt2(sqrt_wt.begin(), l1, 1, false); 
  
  //  NumericVector invwt=1/sqrt(wt);
  
  // Moving Loop inside the function is key for speed
  
  NumericVector yy(l1);
  NumericVector res(m1);
  arma::colvec res2(res.begin(),m1,false); // Reuse memory - update both below
  
  for(int i=0;i<m1;i++){
    
    // Grab one column at a time from b and one row at a time from res
    
    b2temp=b(Range(0,l2-1),Range(i,i));
    
    // Point b2 to memory for that column
    
    arma::mat b2(b2temp.begin(), l2, 1, false); 
    arma::mat restemp(res.begin()+i, 1, 1, false); 
    
    // calculate weighted residuals (element by element multiplication with weights)
    
    xb2=(y2-alpha2- x2 * b2)%sqrt_wt2;
    
    // This is where RSS should be calculated
    // Not sure if this will complain about type differences
    
    restemp=trans(xb2)*xb2;
    
  }
  
  return res;      
  
}








// [[Rcpp::export("rss_face_at_disp")]]

double rss_face_at_disp(double dispersion,
                               Rcpp::List cache,
                               Rcpp::NumericVector cbars_j,
                               Rcpp::NumericVector y,
                               Rcpp::NumericMatrix x,
                               Rcpp::NumericVector alpha,
                               Rcpp::NumericVector wt) {
  // Build 1×l1 matrix, then transpose to l1×1 for Inv_f3_with_disp
  int l1 = cbars_j.size();
  Rcpp::NumericMatrix cbars_small(1, l1);
  for (int k = 0; k < l1; ++k) cbars_small(0, k) = cbars_j[k];
  
  arma::mat theta_row = Inv_f3_with_disp(cache, dispersion, Rcpp::transpose(cbars_small));
  arma::vec beta = theta_row.t(); // 1×l1 -> l1×1
  
  arma::vec y2(y.begin(), y.size(), false);
  arma::vec a2(alpha.begin(), alpha.size(), false);
  arma::mat X(x.begin(), x.nrow(), x.ncol(), false);
  arma::vec w(wt.begin(), wt.size(), false);
  
  arma::vec resid = (y2 - a2 - X * beta) % arma::sqrt(w);
  return arma::as_scalar(resid.t() * resid);
}



// [[Rcpp::export]]
double UB2(double dispersion,
           Rcpp::List cache,
           Rcpp::NumericVector cbars_j,
           Rcpp::NumericVector y,
           Rcpp::NumericMatrix x,
           Rcpp::NumericVector alpha,
           Rcpp::NumericVector wt,
           double rss_min_global) {
  
  // Call the existing RSS function
  double rss_val = rss_face_at_disp(dispersion, cache, cbars_j, y, x, alpha, wt);
  
  // Compute UB2
  double UB2_val = (0.5 / dispersion) * (rss_val - rss_min_global);
  
  return UB2_val;
}


// Utility: safe max for NumericVector
static inline double max_vec(const NumericVector& v) {
  double m = R_NegInf;
  for (int i = 0; i < v.size(); ++i) if (v[i] > m) m = v[i];
  return m;
}


NumericVector EnvBuildLinBound_cpp(NumericMatrix thetabars,
                                   NumericMatrix cbars,
                                   NumericVector y,
                                   NumericMatrix x,
                                   NumericMatrix P,
                                   NumericVector alpha,
                                   double dispstar) {
  // Convert to Armadillo
  arma::mat thetabarsA = as<arma::mat>(thetabars);
  arma::mat cbarsA     = as<arma::mat>(cbars);
  arma::vec yA         = as<arma::vec>(y);
  arma::mat xA         = as<arma::mat>(x);
  arma::mat PA         = as<arma::mat>(P);
  arma::vec alphaA     = as<arma::vec>(alpha);
  
  int gs = cbarsA.n_rows;
  
  arma::mat XtX   = xA.t() * xA;
  arma::vec rhs   = xA.t() * (yA - alphaA);
  arma::mat M     = XtX + dispstar * PA;
  arma::mat Minv  = arma::inv(M);           // match R's solve(M)
  arma::mat H1    = -Minv * PA * Minv;
  
  arma::mat V = -thetabarsA * PA + cbarsA;                 // gs x p
  arma::mat Minv_cbars = cbarsA * Minv.t();                // gs x p
  arma::vec term1 = arma::sum(V % Minv_cbars, 1);
  
  arma::mat rhs_mat = arma::repmat(rhs, 1, gs);            // p x gs
  arma::mat H1_rhs  = (H1 * (rhs_mat + dispstar * cbarsA.t())).t(); // gs x p
  arma::vec term2 = arma::sum(V % H1_rhs, 1);
  
  arma::vec result = term1 + term2;
  
  // Return explicitly as NumericVector
  NumericVector out(gs);
  std::copy(result.begin(), result.end(), out.begin());
  return out;
}


NumericVector thetabar_const_cpp(NumericMatrix P,
                                 NumericMatrix cbars,
                                 NumericMatrix thetabars) {
  arma::mat PA         = as<arma::mat>(P);
  arma::mat cbarsA     = as<arma::mat>(cbars);
  arma::mat thetabarsA = as<arma::mat>(thetabars);
  
  int gs = cbarsA.n_rows;
  arma::vec thetaconst(gs);
  
  for (int j = 0; j < gs; ++j) {
    arma::vec theta_temp = thetabarsA.row(j).t();
    arma::vec cbars_temp = cbarsA.row(j).t();
    thetaconst[j] = -0.5 * arma::as_scalar(theta_temp.t() * PA * theta_temp)
      + arma::as_scalar(cbars_temp.t() * theta_temp);
  }
  
  NumericVector out(gs);
  std::copy(thetaconst.begin(), thetaconst.end(), out.begin());
  return out;
}


// --- Internal helper: RSS pilot timing block ---
// Not exported to R
Rcpp::List run_rss_pilot_block(const Rcpp::Function& parallel_fn,
                               int gs, int l1,
                               double low, double upp,
                               const Rcpp::List& cache,
                               const Rcpp::NumericMatrix& cbars,
                               const Rcpp::NumericVector& y,
                               const Rcpp::NumericMatrix& x,
                               const Rcpp::NumericVector& alpha,
                               const Rcpp::NumericVector& wt,
                               bool use_parallel,
                               bool verbose) {
  double est_total = 0.0;
  // const int pilot_threshold = static_cast<int>(std::pow(3, 10)); // 59,049 faces
  
    // --- Warm-up pilot size ---
    int k1 = std::min(gs, 500);
    
    // Fractional pilots: ~0.5% and ~1.0% of total faces
    auto frac_round = [](double v) { return static_cast<int>(std::round(v)); };
    int k2_target = frac_round(0.005 * static_cast<double>(gs));
    int k3_target = frac_round(0.010 * static_cast<double>(gs));
    
    // Floors/caps
    int floor_k2 = 3000, floor_k3 = 6000;
    int cap_k2   = 50000, cap_k3 = 100000;
    
    int k2 = std::min(gs, std::max(floor_k2, std::min(k2_target, cap_k2)));
    int k3 = std::min(gs, std::max(floor_k3, std::min(k3_target, cap_k3)));
    
    if (k2 <= k1) k2 = std::min(gs, std::max(k1 + 1, floor_k2));
    if (k3 <= k2) k3 = std::min(gs, std::max(k2 + 1, floor_k3));
    
    auto make_slice = [&](int k) {
      Rcpp::NumericMatrix cbars_slice(k, l1);
      for (int i = 0; i < k; ++i)
        for (int j = 0; j < l1; ++j)
          cbars_slice(i, j) = cbars(i, j);
      return cbars_slice;
    };
    
    auto now_num = []() {
      return Rcpp::as<double>(
        Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
      );
    };
    
    // Pilot timings
    double t0 = now_num();
    parallel_fn(Rcpp::Named("par0") = 0.5 * (low + upp),
                Rcpp::Named("low") = low,
                Rcpp::Named("upp") = upp,
                Rcpp::Named("cache") = cache,
                Rcpp::Named("cbars") = make_slice(k1),
                Rcpp::Named("y") = y,
                Rcpp::Named("x") = x,
                Rcpp::Named("alpha") = alpha,
                Rcpp::Named("wt") = wt,
                Rcpp::Named("use_parallel") = use_parallel);
    double t1 = now_num();
    double elapsed1 = t1 - t0;
    
    double t2 = now_num();
    parallel_fn(Rcpp::Named("par0") = 0.5 * (low + upp),
                Rcpp::Named("low") = low,
                Rcpp::Named("upp") = upp,
                Rcpp::Named("cache") = cache,
                Rcpp::Named("cbars") = make_slice(k2),
                Rcpp::Named("y") = y,
                Rcpp::Named("x") = x,
                Rcpp::Named("alpha") = alpha,
                Rcpp::Named("wt") = wt,
                Rcpp::Named("use_parallel") = use_parallel);
    double t3 = now_num();
    double elapsed2 = t3 - t2;
    
    double t4 = now_num();
    parallel_fn(Rcpp::Named("par0") = 0.5 * (low + upp),
                Rcpp::Named("low") = low,
                Rcpp::Named("upp") = upp,
                Rcpp::Named("cache") = cache,
                Rcpp::Named("cbars") = make_slice(k3),
                Rcpp::Named("y") = y,
                Rcpp::Named("x") = x,
                Rcpp::Named("alpha") = alpha,
                Rcpp::Named("wt") = wt,
                Rcpp::Named("use_parallel") = use_parallel);
    double t5 = now_num();
    double elapsed3 = t5 - t4;
    
    double denom   = static_cast<double>(k3 - k2);
    double t_face  = (elapsed3 - elapsed2) / std::max(1.0, denom);
    double t_fixed = elapsed1;
    est_total      = t_fixed + static_cast<double>(gs) * t_face;
    
    auto fmt_hms = [](double seconds) {
      int s = static_cast<int>(std::round(seconds));
      int h = s / 3600; s %= 3600;
      int m = s / 60;   s %= 60;
      std::ostringstream oss;
      if (h) oss << h << "h ";
      if (h || m) oss << m << "m ";
      oss << s << "s";
      return oss.str();
    };
    
    Rcpp::Rcout << "[EnvelopeDispersionBuild:RSS:Pilot] k1=" << k1
                << " (" << (100.0 * k1 / gs) << "%) elapsed=" << elapsed1 << "s; "
                << "k2=" << k2 << " (" << (100.0 * k2 / gs) << "%) elapsed=" << elapsed2 << "s; "
                << "k3=" << k3 << " (" << (100.0 * k3 / gs) << "%) elapsed=" << elapsed3 << "s.\n";
    
    Rcpp::Rcout << "[EnvelopeDispersionBuild:RSS:Pilot] t_fixed=" << t_fixed
                << "s, t_face=" << t_face << "s/face.\n";
    
    Rcpp::Rcout << "[EnvelopeDispersionBuild:RSS:Pilot] Estimated full run = "
                << fmt_hms(est_total) << " (" << est_total << "s).\n";
 
  
  return Rcpp::List::create(Rcpp::Named("est_total") = est_total);
}



// --- Internal helper: UB2 pilot timing block ---
// Not exported to R
Rcpp::List run_ub2_pilot_block(const Rcpp::Function& ub2_parallel_fn,
                               int gs, int l1,
                               double low, double upp,
                               const Rcpp::List& cache,
                               const Rcpp::NumericMatrix& cbars,
                               const Rcpp::NumericVector& y,
                               const Rcpp::NumericMatrix& x,
                               const Rcpp::NumericVector& alpha,
                               const Rcpp::NumericVector& wt,
                               double rss_min_global,
                               bool verbose) {
  double est_total = 0.0;
  // const int pilot_threshold = static_cast<int>(std::pow(3, 10)); // 59,049 faces
  
  // --- Warm-up pilot size ---
  int k1 = std::min(gs, 500);
  
  // Fractional pilots: ~0.5% and ~1.0% of total faces
  auto frac_round = [](double v) { return static_cast<int>(std::round(v)); };
  int k2_target = frac_round(0.005 * static_cast<double>(gs));
  int k3_target = frac_round(0.010 * static_cast<double>(gs));
  
  // Floors/caps
  int floor_k2 = 3000, floor_k3 = 6000;
  int cap_k2   = 50000, cap_k3 = 100000;
  
  int k2 = std::min(gs, std::max(floor_k2, std::min(k2_target, cap_k2)));
  int k3 = std::min(gs, std::max(floor_k3, std::min(k3_target, cap_k3)));
  if (k2 <= k1) k2 = std::min(gs, std::max(k1 + 1, floor_k2));
  if (k3 <= k2) k3 = std::min(gs, std::max(k2 + 1, floor_k3));
  
  auto make_slice = [&](int k) {
    Rcpp::NumericMatrix cbars_slice(k, l1);
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < l1; ++j)
        cbars_slice(i, j) = cbars(i, j);
    return cbars_slice;
  };
  
  auto now_num = []() {
    return Rcpp::as<double>(
      Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
    );
  };
  
  // Pilot timings
  double t0 = now_num();
  ub2_parallel_fn(Rcpp::Named("par0")   = 0.5 * (low + upp),
                  Rcpp::Named("low")    = low,
                  Rcpp::Named("upp")    = upp,
                  Rcpp::Named("cache")  = cache,
                  Rcpp::Named("cbars")  = make_slice(k1),
                  Rcpp::Named("y")      = y,
                  Rcpp::Named("x")      = x,
                  Rcpp::Named("alpha")  = alpha,
                  Rcpp::Named("wt")     = wt,
                  Rcpp::Named("rss_min_global") = rss_min_global);
  double t1 = now_num();
  double elapsed1 = t1 - t0;
  
  double t2 = now_num();
  ub2_parallel_fn(Rcpp::Named("par0")   = 0.5 * (low + upp),
                  Rcpp::Named("low")    = low,
                  Rcpp::Named("upp")    = upp,
                  Rcpp::Named("cache")  = cache,
                  Rcpp::Named("cbars")  = make_slice(k2),
                  Rcpp::Named("y")      = y,
                  Rcpp::Named("x")      = x,
                  Rcpp::Named("alpha")  = alpha,
                  Rcpp::Named("wt")     = wt,
                  Rcpp::Named("rss_min_global") = rss_min_global);
  double t3 = now_num();
  double elapsed2 = t3 - t2;
  
  double t4 = now_num();
  ub2_parallel_fn(Rcpp::Named("par0")   = 0.5 * (low + upp),
                  Rcpp::Named("low")    = low,
                  Rcpp::Named("upp")    = upp,
                  Rcpp::Named("cache")  = cache,
                  Rcpp::Named("cbars")  = make_slice(k3),
                  Rcpp::Named("y")      = y,
                  Rcpp::Named("x")      = x,
                  Rcpp::Named("alpha")  = alpha,
                  Rcpp::Named("wt")     = wt,
                  Rcpp::Named("rss_min_global") = rss_min_global);
  double t5 = now_num();
  double elapsed3 = t5 - t4;
  
  // Estimate per-face slope
  double denom   = static_cast<double>(k3 - k2);
  double t_face  = (elapsed3 - elapsed2) / std::max(1.0, denom);
  double t_fixed = elapsed1;
  est_total      = t_fixed + static_cast<double>(gs) * t_face;
  
  auto fmt_hms = [](double seconds) {
    int s = static_cast<int>(std::round(seconds));
    int h = s / 3600; s %= 3600;
    int m = s / 60;   s %= 60;
    std::ostringstream oss;
    if (h) oss << h << "h ";
    if (h || m) oss << m << "m ";
    oss << s << "s";
    return oss.str();
  };
  
  Rcpp::Rcout << "[EnvelopeDispersionBuild:UB2:Pilot] k1=" << k1
              << " (" << (100.0 * k1 / gs) << "%) elapsed=" << elapsed1 << "s; "
              << "k2=" << k2 << " (" << (100.0 * k2 / gs) << "%) elapsed=" << elapsed2 << "s; "
              << "k3=" << k3 << " (" << (100.0 * k3 / gs) << "%) elapsed=" << elapsed3 << "s.\n";
  
  Rcpp::Rcout << "[EnvelopeDispersionBuild:UB2:Pilot] t_fixed=" << t_fixed
              << "s, t_face=" << t_face << "s/face.\n";
  
  Rcpp::Rcout << "[EnvelopeDispersionBuild:UB2:Pilot] Estimated full run = "
              << fmt_hms(est_total) << " (" << est_total << "s).\n";
  
  return Rcpp::List::create(Rcpp::Named("est_total") = est_total);
}



// [[Rcpp::export]]
List EnvelopeDispersionBuild_cpp(
    List Env,
    double Shape,
    double Rate,
    NumericMatrix P,
    NumericVector y,
    NumericMatrix x,
    NumericVector alpha,
    int n_obs,
    double RSS_post,
    double RSS_ML,
    NumericMatrix mu,         // ← new
    NumericVector wt,         // ← new
    double max_disp_perc ,
    Nullable<double> disp_lower ,
    Nullable<double> disp_upper ,
    bool verbose ,
    bool use_parallel    // ← add flag here
  
)
  {
  
  
  // --- NEW: selector for RSS source ---
  // 1 = use minimization (default)
  // 2 = use RSS_ML (skip minimization)
  int RSS_Min_Type = 1;  // change manually for testing
  int UB2_Min_Type = 1;  // change manually for testing
  
  
  // Step 1: Posterior Gamma parameters (precision prior)
  double shape2 = Shape + static_cast<double>(n_obs) / 2.0;
  double rate3  = Rate  + RSS_post / 2.0;
  
  // Step 2: Dispersion bounds (on sigma^2)
  double low, upp;
  if (disp_lower.isNull() || disp_upper.isNull()) {
    // Call R's qgamma for tail quantiles, then invert to get sigma^2 bounds
    Function qgamma("qgamma");
    NumericVector q_low = qgamma(
      Named("p")     = max_disp_perc,
      Named("shape") = shape2,
      Named("rate")  = rate3
    );
    NumericVector q_upp = qgamma(
      Named("p")     = 1.0 - max_disp_perc,
      Named("shape") = shape2,
      Named("rate")  = rate3
    );
    low = 1.0 / q_low[0];
    upp = 1.0 / q_upp[0];
  } else {
    low = as<double>(disp_lower);
    upp = as<double>(disp_upper);
    if (!R_finite(low) || !R_finite(upp))
      stop("disp_lower/disp_upper must be finite.");
    if (low <= 0.0 || upp <= 0.0)
      stop("disp_lower/disp_upper must be positive.");
    if (upp <= low)
      stop("disp_upper must be strictly greater than disp_lower.");
  }
  
  // Step 3: Extract envelope faces
  NumericMatrix cbars     = Env["cbars"];      // gs x l1
  NumericMatrix thetabars = Env["thetabars"];  // gs x l1 (grid of tangencies)
  NumericVector logP1     = Env["logP"];       // length gs
  int gs = cbars.nrow();
  int l1 = cbars.ncol();
  
  /// Step 3B: Precompute elements for finding inverse function for cbars
  
  
  Rcpp::List cache = Inv_f3_precompute_disp(cbars, y, x, mu, P, alpha, wt);
  
  // Step 3C: Minimize RSS over dispersion for each face (optional diagnostics / UB2 prep)
  // Strategy A (pure C++): call a Brent/golden-section minimizer using rss_face_at_disp()
  // Strategy B (R-side): call optim("Brent") on [low, upp] — easier to prototype
  
  // Step 3C: Minimize RSS over dispersion for each face
  Rcpp::Function optim("optim");
  
  // --- NEW: Call parallel helper and time it ---
//  Rcpp::Environment ns3 = Rcpp::Environment::namespace_env("glmbayes");
//  Rcpp::Function rss_fn = ns3["rss_face_at_disp"];
  
//  Rcpp::Function grad_fn("drss_ddisp");   // exported gradient
  
  
  
  // Optionally: keep the best across faces
   double rss_min_global = R_PosInf;
   [[maybe_unused]] double disp_min_global = NA_REAL;
   [[maybe_unused]] int j_best = -1;
  // Extract parallel results
  Rcpp::NumericVector disp_min_parallel(gs);
  Rcpp::NumericVector rss_min_parallel(gs); 
  
    
  if(RSS_Min_Type==1){
  
  if (verbose) {
    // Print total number of faces before entering the loop
    Rcpp::Rcout << "[EnvelopeDispersionBuild] Total faces to process: "
                << gs << "\n";
    
    Rcpp::Function fmt("format");
    Rcpp::Function systime("Sys.time");
    Rcpp::CharacterVector now = fmt(systime(), Rcpp::Named("format") = "%H:%M:%S");
    Rcpp::Rcout << "[EnvelopeDispersionBuild] >>> Starting RSS minimization loop at "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
    
    
  }

  
  
  
    
  // --- NEW: Call parallel helper and time it ---
//  Rcpp::Function parallel_fn("EnvelopeDispersionBuild_parallel");

  // --- NEW: Call parallel helper and time it ---
  Rcpp::Environment ns = Rcpp::Environment::namespace_env("glmbayes");
  Rcpp::Function parallel_fn = ns["EnvelopeDispersionBuild_parallel"];
  
  double est_total = 0.0;  // declare before pilot block
  
    

  // --- Threshold for pilot runs ---
  const int pilot_threshold = static_cast<int>(std::pow(3, 14)); // 59,049 faces


  // --- Conditional run of pilot block ---
  if (gs >= pilot_threshold) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild] Running RSS pilot block (faces="
                << gs << " >= threshold=" << pilot_threshold << ").\n";

    Rcpp::List pilot_res = run_rss_pilot_block(parallel_fn, gs, l1,
                                               low, upp, cache, cbars,
                                               y, x, alpha, wt,
                                               use_parallel, verbose);
    est_total = pilot_res["est_total"];

    if (verbose) {
      Rcpp::Rcout << "[EnvelopeDispersionBuild] run_rss_pilot_block completed; "
                  << "est_total=" << est_total << " seconds.\n";
    }
  }
   else {
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeDispersionBuild] Skipping RSS pilot block "
                  << "(faces=" << gs << " < threshold=" << pilot_threshold << ").\n";
    }
         

   }  
  
    // --- After computing est_total ---
     double est_total_sec = est_total;  // from pilot estimate
  
  // --- yes/no option if estimate exceeds 5 minutes ---
  if (est_total_sec > 300.0) {
    std::string prompt = "Estimated minimization exceeds 5 minutes. Continue? [y/N]: ";
    
    Rcpp::Function r_interactive("interactive");
    bool is_interactive = Rcpp::as<bool>(r_interactive());
    
    if (is_interactive) {
      Rcpp::Function readline("readline");
      while (true) {
        std::string ans = Rcpp::as<std::string>(readline(Rcpp::wrap(prompt)));
        // trim whitespace
        auto ltrim = [](std::string &s) {
          s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                          [](unsigned char ch){ return !std::isspace(ch); }));
        };
        auto rtrim = [](std::string &s) {
          s.erase(std::find_if(s.rbegin(), s.rend(),
                               [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
        };
        ltrim(ans); rtrim(ans);
        
        if (ans == "y" || ans == "yes" || ans == "1" || ans == "continue") {
          Rcpp::Rcout << "[INFO] User chose to continue full run.\n";
          Rcpp::Rcout << ">>> Running Full parallel Minimization: "
                      << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
                      << "\n";
          break; // proceed to parallel Minimization
        } else if (ans == "n" || ans == "no" || ans == "2" || ans.empty()) {
          Rcpp::Rcout << "[INFO] User declined. Stopping Minimization.\n";
          Rcpp::stop("Minimization stopped by user after time estimate.");
        } else {
          Rcpp::Rcout << "Invalid input. Please enter y (continue) or N (stop).\n";
        }
      }
    } else {
      // Non-interactive session (e.g. CI/CRAN): auto-approve
      Rcpp::Rcout << "[NOTE] Non-interactive session: proceeding automatically.\n";
      Rcpp::Rcout << "[INFO] Proceeding with full run.\n";
      Rcpp::Rcout << ">>> Running Full parallel Minimization: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
                  << "\n";
    }
  }  
    
    
    
  double start_time_parallel = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  
  Rcpp::List parallel_res = parallel_fn(
    Rcpp::Named("par0")   = 0.5 * (low + upp),
    Rcpp::Named("low")    = low,
    Rcpp::Named("upp")    = upp,
    Rcpp::Named("cache")  = cache,
    Rcpp::Named("cbars")  = cbars,
    Rcpp::Named("y")      = y,
    Rcpp::Named("x")      = x,
    Rcpp::Named("alpha")  = alpha,
    Rcpp::Named("wt")     = wt
  );
  
  double end_time_parallel = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  double elapsed_parallel = end_time_parallel - start_time_parallel;
  
  // Break elapsed into h/m/s
  int h_elapsed = static_cast<int>(elapsed_parallel / 3600);
  int m_elapsed = static_cast<int>((elapsed_parallel - h_elapsed*3600) / 60);
  int s_elapsed = static_cast<int>(elapsed_parallel - h_elapsed*3600 - m_elapsed*60);
  
  
  
    // Extract parallel results
     disp_min_parallel = parallel_res["disp_min"];
     rss_min_parallel  = parallel_res["rss_min"];
  

  if (verbose) {
    Rcpp::Function fmt("format");
    Rcpp::Function systime("Sys.time");
    Rcpp::CharacterVector now = fmt(systime(), Rcpp::Named("format") = "%H:%M:%S");
    Rcpp::Rcout << "[EnvelopeDispersionBuild] >>> Exiting RSS minimization loop at "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
    Rcpp::Rcout << "[EnvelopeDispersionBuild] RSS Parallel helper completed in "
                << h_elapsed << "h " << m_elapsed << "m " << s_elapsed << "s.\n";  
    
    
      }

  for (int j = 0; j < gs; ++j) {
    if (rss_min_parallel[j] < rss_min_global) {
      rss_min_global = rss_min_parallel[j];
      disp_min_global = disp_min_parallel[j];
      j_best = j;
    }
  }  
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild] RSS source = MIN (optimized)\n";
  }
  
  
  }
  else { // RSS_Min_Type == 2
    rss_min_global = RSS_ML;
    disp_min_global = 0.5*(low+upp);
    j_best = -1;
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeDispersionBuild] RSS source = ML (skip minimization)\n";
      Rcpp::Rcout << "[EnvelopeDispersionBuild] RSS_ML = " << RSS_ML << "\n";
    }
  }
  
  
  
  
///   (1/low)*(rss_min_parallel[j]-rss_min_global))


  //////////////////////////////////////

  if (verbose) {
    Rcpp::Function fmt("format");
    Rcpp::Function systime("Sys.time");
    Rcpp::CharacterVector now = fmt(systime(), Rcpp::Named("format") = "%H:%M:%S");
    Rcpp::Rcout << "[EnvelopeDispersionBuild] >>> Starting UB2 minimization loop at "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
    
    
  }
  
  
  /// Switch to using namesspace
  
//  Rcpp::Environment ns = Rcpp::Environment::namespace_env("glmbayes");
//  Rcpp::Function ub2_fn = ns["UB2"];
  
  // Preallocate to gs faces
  Rcpp::NumericVector disp_min_ub2(gs);
  Rcpp::NumericVector ub2_min(gs);
  
  if(UB2_Min_Type==1){
    
  
  // --- NEW: Call UB2 parallel helper and time it ---
  Rcpp::Environment ns2 = Rcpp::Environment::namespace_env("glmbayes");
    
  Rcpp::Function ub2_parallel_fn = ns2["EnvelopeUB2_parallel"];
  
  double est_total_ub2 = 0.0;  // declare before pilot block
  
  // --- Threshold for UB2 pilot runs ---
  const int pilot_threshold_ub2 = static_cast<int>(std::pow(3, 14)); // 4,782,969 faces
  
  // --- Conditional run of UB2 pilot block ---
  if (gs >= pilot_threshold_ub2) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild] Running UB2 pilot block (faces="
                << gs << " >= threshold=" << pilot_threshold_ub2 << ").\n";
    
    Rcpp::List ub2_res = run_ub2_pilot_block(ub2_parallel_fn, gs, l1,
                                             low, upp, cache, cbars,
                                             y, x, alpha, wt,
                                             rss_min_global,
                                             verbose);
    est_total_ub2 = ub2_res["est_total"];
    
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeDispersionBuild] run_ub2_pilot_block completed; "
                  << "est_total=" << est_total_ub2 << " seconds.\n";
    }
  } else {
    if (verbose) {
      Rcpp::Rcout << "[EnvelopeDispersionBuild] Skipping UB2 pilot block "
                  << "(faces=" << gs << " < threshold=" << pilot_threshold_ub2 << ").\n";
    }
  }
  
    if (est_total_ub2 > 300.0) {
    std::string prompt = "Estimated UB2 minimization exceeds 5 minutes. Continue? [y/N]: ";
    Rcpp::Function r_interactive("interactive");
    bool is_interactive = Rcpp::as<bool>(r_interactive());
    
    if (is_interactive) {
      Rcpp::Function readline("readline");
      while (true) {
        std::string ans = Rcpp::as<std::string>(readline(Rcpp::wrap(prompt)));
        // trim whitespace
        ans.erase(ans.begin(), std::find_if(ans.begin(), ans.end(),
                            [](unsigned char ch){ return !std::isspace(ch); }));
        ans.erase(std::find_if(ans.rbegin(), ans.rend(),
                               [](unsigned char ch){ return !std::isspace(ch); }).base(), ans.end());
        
        if (ans == "y" || ans == "yes" || ans == "1" || ans == "continue") {
          Rcpp::Rcout << "[INFO] User chose to continue UB2 minimization.\n";
          Rcpp::Rcout << ">>> Running Full UB2 parallel minimization: "
                      << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
                      << "\n";
          break;
        } else if (ans == "n" || ans == "no" || ans == "2" || ans.empty()) {
          Rcpp::Rcout << "[INFO] User declined. Stopping UB2 minimization.\n";
          Rcpp::stop("UB2 minimization stopped by user after time estimate.");
        } else {
          Rcpp::Rcout << "Invalid input. Please enter y (continue) or N (stop).\n";
        }
      }
    } else {
      Rcpp::Rcout << "[NOTE] Non-interactive session: proceeding automatically.\n";
      Rcpp::Rcout << "[INFO] Proceeding with full UB2 minimization.\n";
      Rcpp::Rcout << ">>> Running Full UB2 parallel minimization: "
                  << Rcpp::as<std::string>(Rcpp::Function("format")(Rcpp::Function("Sys.time")()))
                  << "\n";
    }
  }
  
  // --- Run full UB2 minimization ---
  double start_time_ub2 = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  
  // --- Print rss_min_global before UB2 minimization call ---
  if (verbose) {
    Rcpp::Rcout << "[EnvelopeDispersionBuild] rss_min_global_used in optimization is: "
                << rss_min_global << "\n";
  }
  
  
  Rcpp::List ub2_parallel_res = ub2_parallel_fn(
    Rcpp::Named("par0")   = 0.5 * (low + upp),
    Rcpp::Named("low")    = low,
    Rcpp::Named("upp")    = upp,
    Rcpp::Named("cache")  = cache,
    Rcpp::Named("cbars")  = cbars,
    Rcpp::Named("y")      = y,
    Rcpp::Named("x")      = x,
    Rcpp::Named("alpha")  = alpha,
    Rcpp::Named("wt")     = wt,
    Rcpp::Named("rss_min_global") = rss_min_global
  );
  
  double end_time_ub2 = Rcpp::as<double>(
    Rcpp::Function("as.numeric")(Rcpp::Function("Sys.time")())
  );
  
  double elapsed_ub2 = end_time_ub2 - start_time_ub2;
  
  // Break elapsed into h/m/s
  int h_elapsed_ub2 = static_cast<int>(elapsed_ub2 / 3600);
  int m_elapsed_ub2 = static_cast<int>((elapsed_ub2 - h_elapsed_ub2*3600) / 60);
  int s_elapsed_ub2 = static_cast<int>(elapsed_ub2 - h_elapsed_ub2*3600 - m_elapsed_ub2*60);
  
  
  // Extract UB2 parallel results
   disp_min_ub2 = ub2_parallel_res["disp_min"];
   ub2_min      = ub2_parallel_res["ub2_min"];

   for (int j = 0; j < gs; ++j) {
     
   NumericVector cbars_j = cbars(j, _);
  
   }
  
  if (verbose) {
    Rcpp::Function fmt("format");
    Rcpp::Function systime("Sys.time");
    Rcpp::CharacterVector now = fmt(systime(), Rcpp::Named("format") = "%H:%M:%S");
    Rcpp::Rcout << "[EnvelopeDispersionBuild] >>> Exiting UB2 minimization loop at "
                << Rcpp::as<std::string>(now[0]) << " <<<\n";
    Rcpp::Rcout << "[EnvelopeDispersionBuild] UB2 parallel helper completed in "
                << h_elapsed_ub2 << "h " << m_elapsed_ub2 << "m " << s_elapsed_ub2 << "s.\n";
    
      }
  

  // Find global UB2 minimum
  [[maybe_unused]] double ub2_min_global = R_PosInf;
  [[maybe_unused]] double disp_min_global_ub2 = NA_REAL;
  [[maybe_unused]]  int j_best_ub2 = -1;

  
  for (int j = 0; j < gs; ++j) {
//    Rcpp::Rcout << "Index j: " << j << ", ub2_min[j]: " << ub2_min[j] << ", disp_min_ub2[j]: " << disp_min_ub2[j] << std::endl;
    if (ub2_min[j] < ub2_min_global) {
      ub2_min_global = ub2_min[j];
      disp_min_global_ub2 = disp_min_ub2[j];
      j_best_ub2 = j;
    }
  }
  
  

  }
  else { // UB2_Min_Type == 2
    if (RSS_Min_Type == 1) {
      // RSS minimized, UB2 skipped: derive ub2_min from rss_min_parallel
      for (int j = 0; j < gs; ++j) {
        ub2_min[j]      = (0.5 / upp) * (rss_min_parallel[j] - rss_min_global);
        disp_min_ub2[j] = upp;  // enforce upper bound anchor
      }
      if (verbose) {
        Rcpp::Rcout << "[EnvelopeDispersionBuild] UB2 source = derived from RSS_min (skip UB2)\n";
      }
      
    } else if (RSS_Min_Type == 2) {
      // RSS minimized, UB2 skipped: derive ub2_min from rss_min_parallel
      for (int j = 0; j < gs; ++j) {
        ub2_min[j]      = 0;
        disp_min_ub2[j] = upp;  // enforce upper bound anchor
      }
      if (verbose) {
        Rcpp::Rcout << "[EnvelopeDispersionBuild] UB2 source = Set to 0 (skip RSS_Min and UB2 Min)\n";
      }
      
    }
  }
  
  
  
  // Step 4: Base face constants via R helper (keep R version for now)
//  Function thetabar_const_R("thetabar_const");
//  NumericVector thetabar_const_base =
//    thetabar_const_R(P, cbars, thetabars);   // length gs
  
  NumericVector thetabar_const_base =
    thetabar_const_cpp(P, cbars, thetabars);
  
    
  // Step 5: initial anchor (posterior mean; optional)
  // Note: consistency with external rate2: here rate3=Rate+RSS_post/2 as in R function
  double dispstar = rate3 / (shape2 - 1.0);
  

  // Step 6: Face slopes at dispstar via R helper
//  Function EnvBuildLinBound_R("EnvBuildLinBound");
//  NumericVector New_LL_Slope =
//    EnvBuildLinBound_R(thetabars, cbars, y, x, P, alpha, dispstar); // length gs
  

  NumericVector New_LL_Slope =
    EnvBuildLinBound_cpp(thetabars, cbars, y, x, P, alpha, dispstar);
  

  // Step 7: Linear extrapolation of face constants to bounds
  NumericVector thetabar_const_upp_apprx(gs), thetabar_const_low_apprx(gs);
  for (int j = 0; j < gs; ++j) {
    thetabar_const_upp_apprx[j] = thetabar_const_base[j] + (upp - dispstar) * New_LL_Slope[j];
    thetabar_const_low_apprx[j] = thetabar_const_base[j] + (low - dispstar) * New_LL_Slope[j];
  }
  
  // Step 8: Global upper line geometry (match original mean-slope correction)
  double max_low = max_vec(thetabar_const_low_apprx);
  double max_upp = max_vec(thetabar_const_upp_apprx);
  
  // No-op in original; keep for parity via mean slope correction
  double m_New_LL_Slope = Rcpp::mean(New_LL_Slope);
  double max_low_mean   = max_upp - m_New_LL_Slope * (upp - low);
  max_low = max_low_mean;
  
  double new_slope = (max_upp - max_low) / (upp - low);
  double new_int   = max_low - new_slope * low;
  
  // Step 9a: Dispersion anchor (exactly as in original: b1/(-c1))
  double b1 = (upp - low);
  double c1 = -std::log(upp / low);
  dispstar  = b1 / (-c1);  // equivalently (upp - low)/log(upp/low)
  

  // Step 9: Mixture weights per face (match original)
  NumericVector New_logP2(gs);
  NumericVector prob_factor(gs);
  NumericVector prob_factor2(gs);
  for (int j = 0; j < gs; ++j) {
    
    Rcpp::checkUserInterrupt();  // allow user to break out
    
    // cbars_temp is row j (length l1)
    double norm2 = 0.0;
    for (int k = 0; k < l1; ++k) {
      double cjk = cbars(j, k);
      norm2 += cjk * cjk;
    }
    New_logP2[j] = logP1[j] + 0.5 * norm2;
    
    double pf_upp = thetabar_const_upp_apprx[j] - max_upp;
    double pf_low = thetabar_const_low_apprx[j] - max_low;
    prob_factor[j] = (pf_upp > pf_low ? pf_upp : pf_low);
    prob_factor2[j] =prob_factor[j]-ub2_min[j];

    
  }
  

  // Log-space prob factors (kept separate for UB_list, as in R)
  NumericVector lg_prob_factor = clone(prob_factor);
  NumericVector lg_prob_factor2 = clone(prob_factor2);
  

  // Normalize weights (PLSD)
  NumericVector prob_factor_exp(gs);
  NumericVector prob_factor_exp2(gs);
  for (int j = 0; j < gs; ++j){
    
    Rcpp::checkUserInterrupt();  // allow user to break out
    
    
    
    prob_factor_exp[j] = std::exp(New_logP2[j] + prob_factor[j]);
    prob_factor_exp2[j] = std::exp(New_logP2[j] + prob_factor2[j]);
    
  }
  double sumP = std::accumulate(prob_factor_exp.begin(), prob_factor_exp.end(), 0.0);
  double sumP2 = std::accumulate(prob_factor_exp2.begin(), prob_factor_exp2.end(), 0.0);
  for (int j = 0; j < gs; ++j){
    prob_factor_exp[j] /= sumP;
    prob_factor_exp2[j] /= sumP2;
    
  }   
  // Step 10: Envelope constants for dispersion and gamma tilt
  double lm_log2 = new_slope * dispstar;
  double lm_log1 = new_int + new_slope * dispstar - new_slope * std::log(dispstar);
  double shape3  = shape2 - lm_log2;
  
  // Step 11: Package outputs
//  Env["PLSD"] = prob_factor_exp;
  Env["PLSD"] = prob_factor_exp2;
  
  List gamma_list = List::create(
    Named("shape3")     = shape3,
//    Named("rate2")      = Rate + RSS_ML / 2.0,  // matches original definition
    Named("rate2")      = Rate + rss_min_global / 2.0,  // matches original definition
    Named("disp_upper") = upp,
    Named("disp_lower") = low
  );
  
  List UB_list = List::create(
    Named("RSS_ML")         = RSS_ML,               // not RSS_post
    Named("RSS_Min")        = rss_min_global,       // Minimum across faces
    Named("max_New_LL_UB")  = max_upp,
    Named("max_LL_log_disp")= lm_log1 + lm_log2 * std::log(upp),
    Named("lm_log1")        = lm_log1,
    Named("lm_log2")        = lm_log2,
    Named("lg_prob_factor") = lg_prob_factor,
    Named("lmc1")           = new_int,
    Named("lmc2")           = new_slope,
    Named("UB2min")           = ub2_min
  
  );
  
  List diagnostics = List::create(
    Named("dispstar")     = dispstar,
    Named("New_LL_Slope") = New_LL_Slope,
    Named("shape2")       = shape2,
    Named("rate3")        = rate3,
    Named("shape3")       = shape3,
    Named("max_low")      = max_low,
    Named("max_upp")      = max_upp,
    Named("new_slope")    = new_slope,
    Named("new_int")      = new_int,
    Named("prob_factor")  = prob_factor_exp,
    Named("UB2min")           = ub2_min
//  Named("prob_factor2")  = prob_factor_exp2
  );
  
  if (verbose) {
    Rcout << "EnvelopeDispersionBuild diagnostics:\n";
    Rcout << "  dispstar      = " << dispstar << "\n";
    Rcout << "  new_slope     = " << new_slope << "\n";
    Rcout << "  new_int       = " << new_int << "\n";
    Rcout << "  lm_log1       = " << lm_log1 << "\n";
    Rcout << "  lm_log2       = " << lm_log2 << "\n";
    Rcout << "  shape3        = " << shape3 << "\n";
    Rcout << "  max_low       = " << max_low << "\n";
    Rcout << "  max_upp       = " << max_upp << "\n";
    Rcout << "  RSS_ML       = " << RSS_ML << "\n";
    Rcout << "  RSS_Min       = " << rss_min_global << "\n";
    Rcout << "  disp_lower       = " << low << "\n";
    Rcout << "  disp_upper       = " << upp << "\n";
  }
  
  return List::create(
    Named("Env_out")    = Env,
    Named("gamma_list") = gamma_list,
    Named("UB_list")    = UB_list,
    Named("diagnostics")= diagnostics
  );
}






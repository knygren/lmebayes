#include "RcppArmadillo.h"
using namespace Rcpp;


namespace glmbayes {
namespace sim {

Rcpp::List rNormalGammaReg(
    int n,
    Rcpp::NumericVector y,
    Rcpp::NumericMatrix x,
    Rcpp::NumericVector mu,
    Rcpp::NumericMatrix P,
    Rcpp::NumericVector offset,
    Rcpp::NumericVector wt,
    double shape,
    double rate,
    Rcpp::Nullable<double> max_disp_perc,
    Rcpp::Nullable<double> disp_lower,
    Rcpp::Nullable<double> disp_upper,
    bool verbose
) {
  // --------------------------------------------------------------
  // Shell only — no computation yet
  // --------------------------------------------------------------
  
  int nvars = x.ncol();
  
  // Dummy coefficient matrix (n × nvars)
  Rcpp::NumericMatrix coef(n, nvars);
  
  // Dummy dispersion vector
  Rcpp::NumericVector dispersion(n, 1.0);
  
  // Dummy draws vector
  Rcpp::IntegerVector draws(n, 1);
  
  // Empty envelope list (consistent with other functions)
  Rcpp::List Envelope = Rcpp::List::create();
  
  return Rcpp::List::create(
    Rcpp::Named("coefficients")   = coef,
    Rcpp::Named("coef.mode")      = Rcpp::NumericVector(nvars, 0.0),
    Rcpp::Named("dispersion")     = dispersion,
    Rcpp::Named("offset")         = offset,
    Rcpp::Named("Prior")          = Rcpp::List::create(
      Rcpp::Named("mu")        = mu,
      Rcpp::Named("Precision") = P
    ),
    Rcpp::Named("prior.weights")  = wt,
    Rcpp::Named("y")              = y,
    Rcpp::Named("x")              = x,
    Rcpp::Named("iters")          = draws,
    Rcpp::Named("Envelope")       = Envelope
  );
}

} // namespace sim
} // namespace glmbayes
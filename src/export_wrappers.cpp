#include "RcppArmadillo.h"
#include "openclPort.h"

using namespace openclPort;


// =============================================================================
// Only one R-callable export remains: has_opencl, called via
// .has_opencl_cpp() in R/rcpp_wrappers.R -> has_opencl() in R/has_opencl.R.
//
// All Tier 1-4 wrappers (rNormalGLM, Envelope*, rIndepNormalGammaReg_std,
// glmb_Standardize_Model, etc.) and the remaining Tier 5 OpenCL diagnostics
// (load_kernel_source_wrapper, load_kernel_library_wrapper,
// get_opencl_core_count, gpu_names) were removed: no R callers.
// The underlying C++ functions they wrapped are still present in their
// respective *.cpp translation units and remain callable from C++.
// =============================================================================

// [[Rcpp::export]]
bool has_opencl_cpp_export() {
  return has_opencl();
}


// =============================================================================
// Phased Out (no R wrappers; C++ exports commented out)
// - rss_face_at_disp, UB2: former RSS/UB2 minimization; active path uses
//   closed-form C++ bounds.
//
// To fully remove: delete this block, then (1) remove *.o from src/,
// (2) uninstall old glmbayes, (3) Rcpp::compileAttributes(),
// (4) devtools::document(), (5) devtools::install().
// =============================================================================

/*
// [[Rcpp::export]]
double rss_face_at_disp_cpp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt
) {
  return rss_face_at_disp(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt
  );
}

// [[Rcpp::export]]
double UB2_cpp_export(
    double dispersion,
    const Rcpp::List& cache,
    const Rcpp::NumericVector& cbars_j,
    const Rcpp::NumericVector& y,
    const Rcpp::NumericMatrix& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& wt,
    double rss_min_global
) {
  return UB2(
    dispersion,
    cache,
    cbars_j,
    y,
    x,
    alpha,
    wt,
    rss_min_global
  );
}
*/

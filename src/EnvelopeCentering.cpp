// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include "Envelopefuncs.h"
#include <math.h>

using namespace Rcpp;

namespace {

/// Closed-form E[ sum_i w_i (y_i - x_i' beta - offset_i)^2 ] under the same
/// Gaussian posterior as rNormalReg (weighted Normal–Normal conjugate update).
double RSS_helper(
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double dispersion
) {
  Function asMat("as.matrix");
  const int l1 = x.ncol();
  const int l2 = x.nrow();

  NumericMatrix mu2a = asMat(mu);
  NumericMatrix x2b = clone(x);
  arma::mat x2bb(x2b.begin(), l2, l1, false);
  arma::mat P2(P.begin(), P.nrow(), P.ncol(), false);

  NumericVector wt2 = wt / dispersion;
  NumericVector y1 = y - offset;
  arma::vec y2b(y1.begin(), l2, false);
  NumericMatrix W1(l2 + l1, l1);
  arma::mat W(W1.begin(), W1.nrow(), W1.ncol(), false);
  NumericVector z1(l2 + l1);
  arma::vec z(z1.begin(), l2 + l1, false);

  int i;
  for (i = 0; i < l2; i++) {
    x2b(i, _) = x2b(i, _) * sqrt(wt2[i]);
    y1(i) = y1(i) * sqrt(wt2[i]);
  }

  arma::mat RA = arma::chol(P2);
  W.rows(0, l2 - 1) = x2bb;
  W.rows(l2, l2 + l1 - 1) = RA;
  arma::mat mu2(mu2a.begin(), mu2a.nrow(), mu2a.ncol(), false);
  z.rows(0, l2 - 1) = y2b;
  z.rows(l2, l1 + l2 - 1) = RA * mu2;

  Function lm_fit_fun("lm.fit");
  List fit = lm_fit_fun(_["x"] = W, _["y"] = z);
  NumericMatrix b2a = asMat(fit[0]);

  arma::mat IR = arma::inv(arma::trimatu(arma::chol(arma::trans(W) * W)));
  arma::mat Sigma = IR * arma::trans(IR);
  arma::vec b2(b2a.begin(), static_cast<arma::uword>(b2a.nrow() * b2a.ncol()));

  arma::mat X = as<arma::mat>(x);
  arma::vec Y = as<arma::vec>(y);
  arma::vec off = as<arma::vec>(offset);
  arma::vec wv = as<arma::vec>(wt);

  const arma::vec r = Y - X * b2 - off;
  const double rss_at_mean = arma::dot(wv, r % r);
  const arma::mat XtWX = arma::trans(X) * (arma::diagmat(wv) * X);
  const double trace_term = arma::trace(XtWX * Sigma);
  return rss_at_mean + trace_term;
}

}  // namespace

namespace glmbayes {

namespace env {

List EnvelopeCentering(
    NumericVector y,
    NumericMatrix x,
    NumericVector mu,
    NumericMatrix P,
    NumericVector offset,
    NumericVector wt,
    double shape,
    double rate,
    int Gridtype,
    bool verbose
) {
  (void)verbose;
  (void)Gridtype;
  const int n_rss_iter = 10;
  Rcpp::Function lm_wfit("lm.wfit");

  int n_obs = y.size();
  NumericVector ystar(n_obs);
  for (int i = 0; i < n_obs; i++) {
    ystar[i] = y[i] - offset[i];
  }

  double n_w = 0.0;
  for (int i = 0; i < wt.size(); ++i) n_w += wt[i];

  Rcpp::List fit = lm_wfit(
    Rcpp::_["x"] = x,
    Rcpp::_["y"] = ystar,
    Rcpp::_["w"] = wt
  );

  NumericVector res = fit["residuals"];
  double RSS = 0.0;
  for (int i = 0; i < res.size(); i++) {
    RSS += res[i] * res[i];
  }
  int p = Rcpp::as<int>(fit["rank"]);
  double dispersion2 = RSS / (n_obs - p);

  double RSS_post_expected = NA_REAL;

  for (int j = 0; j < n_rss_iter; ++j) {
    const double RSS_closed = RSS_helper(
        y, x, mu, P, offset, wt, dispersion2
    );

    RSS_post_expected = RSS_closed;

    double shape2 = shape + n_w / 2.0;
    double rate2 = rate + RSS_closed / 2.0;
    dispersion2 = rate2 / (shape2 - 1.0);
  }

  return List::create(
    Named("dispersion") = dispersion2,
    Named("RSS_post") = RSS_post_expected
  );
}

}  // namespace env

}  // namespace glmbayes

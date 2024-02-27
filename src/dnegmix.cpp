// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <dnegmix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
arma::vec dnormix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & mean, 
                  const arma::vec & sd) {
  arma::vec ans = arma::zeros<arma::vec>(x.n_rows);
  for (arma::uword j = 0 ; j < w.n_rows ; j++) {
    ans += w.at(j) * arma::normpdf(x, mean.at(j), sd.at(j));
  }
  return ans;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


double dgammix_1(const double & x,
                 const arma::vec & w,
                 const arma::vec & alpha,
                 const arma::vec & beta) {
  double ans = 0.;
  for (arma::uword k = 0; k < w.n_rows; k++) {
    ans += w.at(k) * exp(log_gampdf_1(x, alpha.at(k), beta.at(k)));
  }
  return ans;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


// [[Rcpp::export]]
arma::vec dgammix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & alpha, 
                  const arma::vec & beta) {
  arma::vec ans = arma::zeros<arma::vec>(x.n_rows);
  for (arma::uword k = 0 ; k < x.n_rows ; k++) {
    ans.at(k) = dgammix_1(x.at(k), w, alpha, beta);
  }
  return ans;
}

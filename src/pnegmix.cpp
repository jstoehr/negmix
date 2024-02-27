// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
arma::vec pnormix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & mean, 
                  const arma::vec & sd) {
  arma::vec ans = arma::zeros<arma::vec>(x.n_rows);
  for (arma::uword j = 0 ; j < w.n_rows ; j++) {
    ans += w.at(j) * arma::normcdf(x, mean.at(j), sd.at(j));
  }
  return ans;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


// [[Rcpp::export]]
arma::vec pgammix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & alpha, 
                  const arma::vec & inv_beta) {
  arma::vec ans = arma::zeros<arma::vec>(x.n_rows);
  for (arma::uword i = 0 ; i < x.n_rows ; i ++) {
    for (arma::uword j = 0 ; j < w.n_rows ; j++) {
      ans.at(i) += w.at(j) * R::pgamma(x.at(i), alpha.at(j), inv_beta.at(j), true, false);
    }
  }
  return ans;
}
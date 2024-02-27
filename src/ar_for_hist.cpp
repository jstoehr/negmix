// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_pair.hpp>
#include <utils.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

double ar_for_hist_1(f_pair & f,
                     const arma::vec & grid,
                     const arma::vec & h,
                     const arma::vec & p_cell,
                     double & counter) {
  
  double cdf_left = f.w_cdf_p(grid.at(0))/f.w_p;
  double cdf_right = f.w_cdf_p(grid.at(grid.n_rows - 1))/f.w_p;
  
  // --- For the latent variable
  arma::uvec index = arma::regspace<arma::uvec>(0, p_cell.n_elem - 1);
  arma::uword z = sample_int_1(index, p_cell);
  
  // --- Containers
  double y, rho = 0.;
  
  if (z == 0) {
    // --- Drawing from left queue
    while (arma::randu() > rho) {
      y = f.q_p(arma::randu() * cdf_left);
      rho = 1. - f.w_pdf_n(y)/f.w_pdf_p(y);
      ++counter;
    }
  } else if (z == p_cell.n_elem - 1) {
    // --- Drawing from right queue
    while (arma::randu() > rho) {
      y = f.q_p(arma::randu(arma::distr_param(cdf_right, 1.)));
      rho = 1. - f.w_pdf_n(y)/f.w_pdf_p(y);
      ++counter;
    } 
  } else {
    while (arma::randu() > rho) {
      y = arma::randu(arma::distr_param(grid.at(z - 1), grid.at(z)));
      rho = f.w_pdf_p(y) - f.w_pdf_n(y);
      rho /= h.at(z - 1);
      ++counter;
    }
  }
  
  return y;
}
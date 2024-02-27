// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void ar_for_model(f_model & f,
                  arma::vec & ans,
                  double & p_accept,
                  double & cpu_time) {
  clock_t start, end;
  start = clock();
  
  // --- For sampling latent variable
  double m = arma::sum(f.w_p);
  arma::vec prob = f.w_p/m;
  arma::uvec index = arma::regspace<arma::uvec>(0, f.w_p.n_elem - 1);
  
  // --- Containers
  arma::uword z, n = ans.n_rows;
  double y, counter = 0.;
  
  for (arma::uword k = 0; k < n; k++) {
    z = sample_int_1(index, prob);
    y = f.rand_1_from_pc(z);
    ++counter;
    while (arma::randu() > 1. - f.w_pdf_n(y)/f.w_pdf_p(y)) {
      z = sample_int_1(index, prob);
      y = f.rand_1_from_pc(z);
      ++counter;
    }
    ans.at(k) = y;
  }
  end = clock();
  cpu_time = (end - start) / static_cast<double>(CLOCKS_PER_SEC);
  p_accept = static_cast<double>(n)/counter;
}
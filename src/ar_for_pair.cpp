// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_pair.hpp>
#include <build_hist.hpp>
#include <ar_for_hist.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

double ar_for_pair_1(f_pair & f,
                     double & counter) {
  double y = f.rand_1();
  ++counter;
  while (arma::randu() > 1. - f.w_pdf_n(y)/f.w_pdf_p(y)) {
    y = f.rand_1();
    ++counter;
  }
  return y;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


double sample_from_pair(f_pair & f,
                        const double & delta, 
                        const double & eps_d,
                        const double & p_accept,
                        arma::uword & is_build,
                        arma::vec & grid,
                        arma::vec & h,
                        arma::vec & p_cell,
                        double & counter,
                        const bool & optim,
                        const bool & use_mono,
                        const double & n_points,
                        const int & maxit, 
                        const double & eps_f, 
                        const double & eps_g) {
  double y;
  if (p_accept < delta) {
    
    if (is_build == 0) {
      // --- Histogram has not been built yet
      // --- Building histogram
      build_hist(f, delta, eps_d, grid, h, p_cell,
                 optim, use_mono, n_points, maxit, eps_f, eps_g);
      is_build = 1;
    }
    
    // --- Sampling from histogram
    y = ar_for_hist_1(f, grid, h, p_cell, counter);
  } else {
    // --- Sampling from positive part
    y = ar_for_pair_1(f, counter);
  }
  return y;
}

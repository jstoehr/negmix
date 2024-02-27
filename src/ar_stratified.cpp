// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void ar_stratified(f_model & f,
                   const double & delta,
                   const double & eps_d,
                   arma::vec & ans,
                   double & p_accept,
                   double & cpu_time,
                   const bool & optim,
                   const bool & use_mono,
                   const double & n_points,
                   const int & maxit, 
                   const double & eps_f, 
                   const double & eps_g,
                   const double & tol_simplex) {
  clock_t start, end;
  start = clock();
  
  // --- Building pairs
  f.build_pairs(delta, tol_simplex);
  
  arma::uword n_pairs = f.w_pair_p.n_rows;
  
  // --- For sampling latent variable
  arma::vec w_p_prop = join_vert(f.c_norm, f.r_p);
  double m_pair = sum(w_p_prop);
  arma::vec prob = w_p_prop/m_pair;
  arma::uvec index = arma::regspace<arma::uvec>(0, w_p_prop.n_elem - 1);
  
  // --- Building histograms
  arma::uvec is_build = arma::zeros<arma::uvec>(n_pairs);
  arma::field<arma::vec> grid(n_pairs);
  arma::field<arma::vec> h(n_pairs);
  arma::field<arma::vec> p_cell(n_pairs);
  
  // --- Containers
  arma::uword z, n = ans.n_rows;
  double y, ratio, counter = 0.;
  
  for (arma::uword k = 0; k < n; k++) {
    ratio = 0.;
    while (arma::randu() > ratio) {
      z = sample_int_1(index, prob);
      
      if (z < n_pairs) {
        y = f.rand_from_pair(z, delta, eps_d, f.p_accept_pair.at(z),
                             is_build.at(z), grid.at(z), h.at(z), p_cell.at(z),
                             counter, optim, use_mono, n_points, 
                             maxit, eps_f, eps_g);
      } else {
        z -= n_pairs;
        y = f.rand_1_from_pc(f.list_r_p.at(z));
        ++counter;
      }
      ratio = f.stratified_accept_ratio(y);
    }
    ans.at(k) = y;
  }
  end = clock();
  cpu_time = (end - start) / static_cast<double>(CLOCKS_PER_SEC);
  p_accept = static_cast<double>(n)/counter;
}
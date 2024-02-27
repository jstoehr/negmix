// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_gammix.hpp>
#include <f_normix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

arma::vec f_model::inv_cdf(const arma::vec & u,
                           const double & precision,
                           const unsigned & n_bins) {
  
  arma::vec q = this->inv_cdf_init(n_bins);
  arma::vec p = this->cdf(q);
  
  // --- Bin containers
  arma::uvec i_temp;
  arma::uword i;
  double p_inf, p_sup;
  double q_inf, q_sup;
  // --- Containers for the search
  double slope;
  double intercep;
  arma::vec q_new(1);
  arma::vec p_new;
  // --- For inserting new elements without search
  arma::vec q_insert, p_insert;
  
  // --- Output
  arma::vec ans(u.n_rows);
  
  // --- Dealing with the left tail
  // ------ We refine if the lowest value of u
  // ------ is lower than the lowest available 
  // ------ value of p
  arma::uword i_min = u.index_min();
  while(p.at(0) - u(i_min) > precision) {
    // --- SHOULD NEVER HAPPEN WHEN STARTING AT 0
    slope = (p.at(1) - p.at(0))/(q.at(1) - q.at(0));
    intercep = p.at(0) - q.at(0) * slope;
    q_new.at(0) = (u(i_min) - intercep)/slope;
    q.insert_rows(0, q_new);
    p.insert_rows(0, this->cdf(q_new));
  }
  
  // --- Dealing with the right tail
  // ------ We refine if the largest value of u
  // ------ is larger than the largest available 
  // ------ value of p
  arma::uword i_max = u.index_max();
  while(u(i_max) - p.back() > precision) {
    slope = arma::as_scalar(arma::diff(p.tail(2))/arma::diff(q.tail(2)));
    intercep = p.back() - q.back() * slope;
    q_new.at(0) = (u(i_max) - intercep)/slope;
    q.insert_rows(q.n_rows, q_new);
    p.insert_rows(p.n_rows, this->cdf(q_new));
  }
  
  for (arma::uword k = 0; k < u.n_rows; k++) {
    i_temp = arma::find(p > u.at(k), 1, "first");
    if (i_temp.is_empty()) {
      // --- Right tail
      ans.at(k) = q.back();
    } else if (i_temp.at(0) == 0) {
      // --- Left tail
      ans.at(k) = q.front();
    } else {
      i = i_temp.at(0);
      p_inf = p.at(i - 1);
      p_sup = p.at(i);
      q_inf = q.at(i - 1);
      q_sup = q.at(i);
      
      // --- Shrinking the interval till we reach precision
      p_insert.empty();
      q_insert.empty();
      while (std::min(u.at(k) - p_inf, p_sup - u.at(k)) > precision) {
        slope = (p_sup - p_inf)/(q_sup - q_inf);
        intercep = p_inf - q_inf * slope;
        q_new.at(0) = (u.at(k) - intercep)/slope;
        p_new = this->cdf(q_new);
        if (u.at(k) < p_new.at(0)) {
          q_sup = q_new.at(0);
          p_sup = p_new.at(0);
        } else {
          q_inf = q_new.at(0);
          p_inf = p_new.at(0);
        }
      }
      
      // --- Computing the output
      if (2.0 * u.at(k) - p_inf < p_sup) {
        ans.at(k) = q_inf;
      } else {
        ans.at(k) = q_sup;
      }
    }
  }
  return ans;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void inv_cdf_method(f_model & f,
                    arma::vec & ans,
                    double & cpu_time,
                    const double & precision,
                    const unsigned & n_bins) {
  clock_t start, end;
  start = clock();
  
  arma::vec u = arma::randu(ans.n_rows);
  ans = f.inv_cdf(u, precision, n_bins);
  
  end = clock();
  cpu_time = (end - start) / static_cast<double>(CLOCKS_PER_SEC);
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


// [[Rcpp::export]]
arma::vec cpp_qnegmix(const arma::vec & x,
                      const Rcpp::List & par,
                      const std::string & family,
                      const double & precision,
                      const unsigned & n_bins) {
  char sig = family[0];
  
  if (sig == 'n') {
    f_normix f(par["w_p"], par["mean_p"], par["sd_p"], 
               par["w_n"], par["mean_n"], par["sd_n"], 
                   0., 1., 1.);
    return f.inv_cdf(x, precision, n_bins);
  } else if (sig == 'g') {
    f_gammix f(par["w_p"], par["shape_p"], par["rate_p"],
               par["w_n"], par["shape_n"], par["rate_n"],
                   0., 1., 1.);
    return f.inv_cdf(x, precision, n_bins);
  } else {
    Rcpp::stop(" Error in rnegmix. Family should be normal or gamma\n");
    std::exit(0);
  }
}
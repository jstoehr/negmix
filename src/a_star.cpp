// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <a_star.hpp>
#include <dnegmix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
double a_star_norm(const double & mean_p, 
                   const double & sd_p, 
                   const double & mean_n, 
                   const double & sd_n,
                   const bool & do_log) {
  // Return the lower bound on "a" for a mixture of the form (a * f - g)/(a - 1)
  double x_opt = x_star_norm(mean_p, sd_p, mean_n, sd_n);
  double d_p = arma::log_normpdf(x_opt, mean_p, sd_p);
  double d_n = arma::log_normpdf(x_opt, mean_n, sd_n);
  if (do_log) {
    return d_n - d_p;
  } else {
    return exp(d_n - d_p);
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


// [[Rcpp::export]]
double a_star_gam(const double & alpha_p, 
                  const double & beta_p, 
                  const double & alpha_n, 
                  const double & beta_n,
                  const bool & do_log) {
  // Return the lower bound on "a" for a mixture of the form (a * f - g)/(a - 1)
  double a;
  if (alpha_p == alpha_n) {
    a = alpha_p * (log(beta_n) - log(beta_p));
  } else {
    double x_opt = (alpha_p - alpha_n)/(beta_p - beta_n);
    double d_p = log_gampdf_1(x_opt, alpha_p, beta_p);
    double d_n = log_gampdf_1(x_opt, alpha_n, beta_n);
    a = d_n - d_p;
  }
  if (do_log) {
    return a;
  } else {
    return exp(a);
  }
}

// [[Rcpp::export]]
double toto(double a, double b) {
  double c = arma::randi<double>(arma::distr_param(a, b));
  return c;
}

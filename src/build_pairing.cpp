// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_gammix.hpp>
#include <f_normix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

// [[Rcpp::export]]
arma::mat cpp_map_pairs(const double & delta, 
                        const Rcpp::List & par,
                        const std::string & family) {
  
  char sig = family[0];
  
  if (sig == 'n') {
    f_normix f(par["w_p"], par["mean_p"], par["sd_p"],
               par["w_n"], par["mean_n"], par["sd_n"],
                   0., 1., 1.);
    
    return f.map_pairs(delta);
    
  } else if (sig == 'g') {
    f_gammix f(par["w_p"], par["shape_p"], par["rate_p"],
               par["w_n"], par["shape_n"], par["rate_n"],
                   0., 1., 1.);
    
    return f.map_pairs(delta);
  } else {
    Rcpp::stop(" Error in rnegmix. Family should be normal or gamma\n");
    std::exit(0);
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


// [[Rcpp::export]]
Rcpp::List build_pairing(const double & delta, 
                         const Rcpp::List & par,
                         const std::string & family,
                         const double & tol_simplex) {
  
  char sig = family[0];
  
  if (sig == 'n') {
    f_normix f(par["w_p"], par["mean_p"], par["sd_p"],
               par["w_n"], par["mean_n"], par["sd_n"],
                   0., 1., 1.);
    
    f.build_pairs(delta, tol_simplex);
    
    return Rcpp::List::create(
      Rcpp::Named("pairs", f.list_pair + 1),
      Rcpp::Named("w_pair_p", f.w_pair_p),
      Rcpp::Named("w_pair_n", f.w_pair_n),
      Rcpp::Named("residuals_p", f.list_r_p + 1),
      Rcpp::Named("r_p", f.r_p),
      Rcpp::Named("residuals_n", f.list_r_n + 1),
      Rcpp::Named("r_n", f.r_n),
      Rcpp::Named("p_accept", f.p_accept_pairing(delta))
    );
  } else if (sig == 'g') {
    f_gammix f(par["w_p"], par["shape_p"], par["rate_p"],
               par["w_n"], par["shape_n"], par["rate_n"],
                   0., 1., 1.);
    
    f.build_pairs(delta, tol_simplex);
    
    return Rcpp::List::create(
      Rcpp::Named("pairs", f.list_pair + 1),
      Rcpp::Named("w_pair_p", f.w_pair_p),
      Rcpp::Named("w_pair_n", f.w_pair_n),
      Rcpp::Named("residuals_p", f.list_r_p + 1),
      Rcpp::Named("r_p", f.r_p),
      Rcpp::Named("residuals_n", f.list_r_n + 1),
      Rcpp::Named("r_n", f.r_n),
      Rcpp::Named("p_accept", f.p_accept_pairing(delta))
    );
  } else {
    Rcpp::stop(" Error in rnegmix. Family should be normal or gamma\n");
    std::exit(0);
  }
}

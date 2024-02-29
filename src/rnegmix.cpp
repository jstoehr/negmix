// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <build_grid.hpp>
#include <build_hist.hpp>
#include <ar_for_model.hpp>
#include <ar_stratified.hpp>
#include <qnegmix.hpp>
#include <f_gammix.hpp>
#include <f_normix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

Rcpp::List get_grid(const double & delta,
                    const std::string & family,
                    const Rcpp::List & par,
                    const Rcpp::List & control) {
  // --- Containers
  arma::vec base;
  arma::ivec mono;
  char sig = family[0];
  
  if (sig == 'n') {
    f_pair_norm f(par["w_p"], par["mean_p"], par["sd_p"],
                  par["w_n"], par["mean_n"], par["sd_n"],
                      0., 1., control["lambda"], 1);
    build_grid(f, delta, control["eps_d"], base, mono,
               control["optim"], control["use_mono"],
                                        control["maxit"], control["eps_f"], control["eps_g"]);
  } else if (sig == 'g') {
    f_pair_gam f(par["w_p"], par["shape_p"], par["rate_p"],
                 par["w_n"], par["shape_n"], par["rate_n"],
                     0., 1., control["lambda"], 1);
    build_grid(f, delta, control["eps_d"], base, mono,
               control["optim"], control["use_mono"],
                                        control["maxit"], control["eps_f"], control["eps_g"]);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("base", base),
    Rcpp::Named("mono", mono)
  );
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


Rcpp::List get_hist(const double & delta,
                    const std::string & family,
                    const Rcpp::List & par,
                    const Rcpp::List & control) {
  // --- Container
  arma::vec grid, h, p_cell;
  char sig = family[0];
  
  if (sig == 'n') {
    f_pair_norm f(par["w_p"], par["mean_p"], par["sd_p"],
                  par["w_n"], par["mean_n"], par["sd_n"],
                      0., 1., control["lambda"], 1);
    build_hist(f, delta, control["eps_d"], grid, h, p_cell, 
               control["optim"], control["use_mono"], control["n_points"], 
                       control["maxit"], control["eps_f"], control["eps_g"]);
  } else if (sig == 'g') {
    f_pair_gam f(par["w_p"], par["shape_p"], par["rate_p"],
                 par["w_n"], par["shape_n"], par["rate_n"],
                     0., 1., control["lambda"], 1);
    build_hist(f, delta, control["eps_d"], grid, h, p_cell, 
               control["optim"], control["use_mono"], control["n_points"], 
                       control["maxit"], control["eps_f"], control["eps_g"]);
  }
  
  return Rcpp::List::create(
    Rcpp::Named("breaks", grid),
    Rcpp::Named("height", h),
    Rcpp::Named("p_bin", p_cell)
  );
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

Rcpp::List rnegmix_model(const unsigned & n,
                         f_model & f,
                         const std::string & method,
                         const Rcpp::List & control) {
  // --- Containers
  arma::vec ans(n);
  double p_accept, cpu_time;
  char sig = method[0];
  
  if (sig == 'v') {
    ar_for_model(f, ans, p_accept, cpu_time);
    
    return Rcpp::List::create(
      Rcpp::Named("sample", ans),
      Rcpp::Named("prob_accept", 1./arma::sum(f.w_p)),
      Rcpp::Named("freq_accept", p_accept),
      Rcpp::Named("time", cpu_time)
    );
  } else if (sig == 's') {
    ar_stratified(f, control["delta"], control["eps_d"],
                  ans, p_accept, cpu_time,
                  control["optim"], control["use_mono"],
                                           control["n_points"], control["maxit"], control["eps_f"], control["eps_g"],
                                                                                                           control["tol_simplex"]);
    
    return Rcpp::List::create(
      Rcpp::Named("sample", ans),
      Rcpp::Named("prob_accept", f.p_accept_pairing(control["delta"])),
      Rcpp::Named("freq_accept", p_accept),
      Rcpp::Named("time", cpu_time),
      Rcpp::Named("pairs", f.list_pair + 1),
      Rcpp::Named("w_pair_p", f.w_pair_p),
      Rcpp::Named("w_pair_n", f.w_pair_n),
      Rcpp::Named("residuals_p", f.list_r_p + 1),
      Rcpp::Named("r_p", f.r_p),
      Rcpp::Named("residuals_n", f.list_r_n + 1),
      Rcpp::Named("r_n", f.r_n)
    );
    
  } else {
    inv_cdf_method(f, ans, cpu_time,
                   control["precision"], control["breaks"]);
    
    return Rcpp::List::create(
      Rcpp::Named("sample", ans),
      Rcpp::Named("time", cpu_time)
    );
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


// [[Rcpp::export]]
Rcpp::List cpp_rnegmix(const unsigned & n,
                       const Rcpp::List & par,
                       const std::string & family,
                       const std::string & method,
                       const Rcpp::List & control) {
  char sig = family[0];

  if (sig == 'n') {
    f_normix f(par["w_p"], par["mean_p"], par["sd_p"],
               par["w_n"], par["mean_n"], par["sd_n"],
                   0., 1., control["lambda"]);
    return rnegmix_model(n, f, method, control);
  } else if (sig == 'g') {
    f_gammix f(par["w_p"], par["shape_p"], par["rate_p"],
               par["w_n"], par["shape_n"], par["rate_n"],
                   0., 1., control["lambda"]);
    return rnegmix_model(n, f, method, control);
  } else {
    Rcpp::stop(" Error in rnegmix. Family should be normal or gamma\n");
    std::exit(0);
  }
}

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
Rcpp::List render_model(arma::uword & k_p,
                        arma::uword & k_n,
                        const std::string & family,
                        const Rcpp::List & par,
                        const bool & by_2,
                        const Rcpp::List & control) {
  char sig = family[0];
  
  if (sig == 'n') {
    f_normix f(k_p, k_n, 0., 1., control["lambda"]);
    
    // --- Init positive weight component parameter
    arma::vec m_p = par["mean_p"];
    arma::vec s_p = par["sd_p"];
    f.mean_p.head(m_p.n_rows) = m_p;
    f.sd_p.head(m_p.n_rows) = s_p;
    
    if (by_2) {
      f.create_benchmark_by_2(
        control["k_p_init"],
               control["k_p_pair"], control["n_single"], control["n_pair"],
                       control["p_min"], control["p_max"], control["p_star"],
                               control["p_common"], control["p_rank"], control["p_limit"],
                                                                              control["maxit_0"], control["eps_0"]);
    } else {
      f.create_benchmark(
        control["k_p_init"],
               control["p_min"], control["p_max"], control["p_star"],
                       control["p_common"], control["p_rank"],
                                                   control["maxit_0"], control["eps_0"],
                                                                              control["maxit"], control["eps_f"], control["eps_g"]);
    }
    
    return Rcpp::List::create(
      Rcpp::Named("w_p", f.w_p),
      Rcpp::Named("mean_p", f.mean_p),
      Rcpp::Named("sd_p", f.sd_p),
      Rcpp::Named("w_n", f.w_n),
      Rcpp::Named("mean_n", f.mean_n),
      Rcpp::Named("sd_n", f.sd_n),
      Rcpp::Named("control", control)
    );
    
  } else if (sig == 'g') {
    f_gammix f(k_p, k_n, 0., 1., control["lambda"]);
    
    // --- Init positive weight component parameter
    arma::vec a_p = par["shape_p"];
    arma::vec b_p = par["rate_p"];
    f.alpha_p.head(a_p.n_rows) = a_p;
    f.beta_p.head(a_p.n_rows) = b_p;
    
    if (by_2) {
      f.create_benchmark_by_2(
        control["k_p_init"],
               control["k_p_pair"], control["n_single"], control["n_pair"],
                       control["p_min"], control["p_max"], control["p_star"],
                               control["p_common"], control["p_rank"], control["p_limit"],
                                                                              control["maxit_0"], control["eps_0"]);
    } else {
      f.create_benchmark(
        control["k_p_init"],
               control["p_min"], control["p_max"], control["p_star"],
                       control["p_common"], control["p_rank"],
                                                   control["maxit_0"], control["eps_0"],
                                                                              control["maxit"], control["eps_f"], control["eps_g"]);
    }
    
    return Rcpp::List::create(
      Rcpp::Named("w_p", f.w_p),
      Rcpp::Named("shape_p", f.alpha_p),
      Rcpp::Named("rate_p", f.beta_p),
      Rcpp::Named("w_n", f.w_n),
      Rcpp::Named("shape_n", f.alpha_n),
      Rcpp::Named("rate_n", f.beta_n),
      Rcpp::Named("control", control)
    );
    
  } else {
    Rcpp::stop(" Error in rnegmix. Family should be normal or gamma\n");
    std::exit(0);
  }
}
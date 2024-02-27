#ifndef INST_INCLUDE_F_GAMMIX_HPP_
#define INST_INCLUDE_F_GAMMIX_HPP_

#include <f_model.hpp>
#include <f_pair_gam.hpp>
#include <ar_for_pair.hpp>
#include <pnegmix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_gammix: public f_model {
public:
  f_gammix(
    const arma::uword & k_p_, 
    const arma::uword & k_n_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_model(k_p_, k_n_, inf_, sup_, l_),
  alpha_p(arma::zeros(k_p_)), beta_p(arma::zeros(k_p_)),
  alpha_n(arma::zeros(k_n_)), beta_n(arma::zeros(k_n_)){};
  
  f_gammix(
    const arma::vec & w_p_, 
    const arma::vec & alpha_p_, 
    const arma::vec & beta_p_,
    const arma::vec & w_n_, 
    const arma::vec & alpha_n_, 
    const arma::vec & beta_n_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_model(w_p_, w_n_, inf_, sup_, l_),
  alpha_p(alpha_p_), beta_p(beta_p_), inv_beta_p(1./beta_p_),
  alpha_n(alpha_n_), beta_n(beta_n_), inv_beta_n(1./beta_n_){};
  
  virtual ~f_gammix(){};
  
  // --- Methods for creating benchmark
  void render_pair_from_pc(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & p, 
      const bool & is_star, 
      const bool & has_common,
      const int & maxit_0, 
      const double & eps_0
  );
  
  void render_pair_from_nc(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & p, 
      const bool & is_star, 
      const bool & has_common,
      const int & maxit_0, 
      const double & eps_0
  );
  
  void set_rand_par_p(const arma::uword & k) {
    alpha_p.at(k) = alpha_n.min() * arma::randu();
    beta_p.at(k) = beta_n.min() * arma::randu();
  };
  
  void add_single_comp(const arma::uword & k_1, const arma::uword & k_2);
  
  arma::vec w_valid_mixt(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & w_p_loc,
      const double & w_n_loc,
      const int & maxit_0, 
      const double & eps_0, 
      const int & maxit, 
      const double & eps_f, 
      const double & eps_g
  );
  
  void set_alt_param() {
    inv_beta_p = 1./beta_p;
    inv_beta_n = 1./beta_n;
  };
  
  // --- Sampling methods
  double rand_1_from_pc(const arma::uword & z) {
    return arma::randg(arma::distr_param(alpha_p.at(z), inv_beta_p.at(z)));
  };
  
  double rand_from_pair(
      const arma::uword & z,
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
      const double & eps_g
  ) {
    f_pair_gam f(
        w_pair_p.at(z)/c_norm.at(z), alpha_pair_p.at(z), beta_pair_p.at(z),
        w_pair_n.at(z)/c_norm.at(z), alpha_pair_n.at(z), beta_pair_n.at(z),
        0., 1., lambda, 1
    );
    return sample_from_pair(f, delta, eps_d, p_accept, is_build,
                            grid, h, p_cell, counter,
                            optim, use_mono, n_points, 
                            maxit, eps_f, eps_g);
  };
  
  // --- Negative mixture CDF
  arma::vec cdf(const arma::vec & x) {
    return pgammix(x, w_p, alpha_p, inv_beta_p) - pgammix(x, w_n, alpha_n, inv_beta_n);
  }
  
  // --- Negative mixture quantile function
  arma::vec inv_cdf_init(const unsigned & n_bins) {
    double temp, max = 0.;
    for (arma::uword k = 0; k < k_p; k++) {
      temp = R::qgamma(.9, alpha_p.at(k), inv_beta_p.at(k), true, false);
      if (temp > max) {
        max = temp;
      }
    } 
    return arma::linspace(0, max, n_bins);
  };
  
  // --- Negative mixture PDF
  double w_pdf_p(const double & x) {
    return dgammix_1(x, w_p, alpha_p, beta_p);
  };
  
  double w_pdf_n(const double & x) {
    return dgammix_1(x, w_n, alpha_n, beta_n);
  };
  
  // --- Methods for pairs
  bool is_valid_pair(const arma::uword & i, 
                     const arma::uword & j) {
    return (alpha_n.at(j) >= alpha_p.at(i) and beta_n.at(j) > beta_p.at(i));
  };
  
  double a_star_pair(const arma::uword & i, 
                     const arma::uword & j) {
    return a_star_gam(alpha_p.at(i), beta_p.at(i), alpha_n.at(j), beta_n.at(j), false);
  };
  
  void set_pairs_param() {
    const arma::uvec r_0 = list_pair.row(0).t();
    const arma::uvec r_1 = list_pair.row(1).t();
    // --- Extracting alpha and beta
    alpha_pair_p = alpha_p.elem(r_0);
    beta_pair_p = beta_p.elem(r_0);
    alpha_pair_n = alpha_n.elem(r_1);
    beta_pair_n = beta_n.elem(r_1);
    
    // --- Cleaning positive and negative remainder
    alpha_r_p = alpha_p.elem(list_r_p);
    beta_r_p = beta_p.elem(list_r_p);
    
    alpha_r_n = alpha_n.elem(list_r_n);
    beta_r_n = beta_n.elem(list_r_n);
  };
  
  // --- Method for stratified sampling
  double stratified_accept_ratio(const double & y) {
    double t_prop = this->w_pdf_p(y) - dgammix_1(y, w_pair_n, alpha_pair_n, beta_pair_n);
    double t_n_r = dgammix_1(y, r_n, alpha_r_n, beta_r_n);
    return 1. - t_n_r / t_prop;
  };
  
  // --- 
  double f_grad(Constvec & x, Refvec grad) {return 0.;};
  
public:
  arma::vec alpha_p, beta_p, inv_beta_p, alpha_n, beta_n, inv_beta_n;
  arma::vec alpha_pair_p, beta_pair_p, alpha_pair_n, beta_pair_n;
  arma::vec alpha_r_p, beta_r_p;
  arma::vec alpha_r_n, beta_r_n;
};

#endif /* INST_INCLUDE_F_GAMMIX_HPP_ */
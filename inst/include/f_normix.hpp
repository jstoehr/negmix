#ifndef INST_INCLUDE_F_NORMIX_HPP_
#define INST_INCLUDE_F_NORMIX_HPP_

#include <f_pair_norm.hpp>
#include <f_obj_norm.hpp>
#include <pnegmix.hpp>
#include <a_star.hpp>
#include <ar_for_pair.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_normix: public f_model {
public:
  f_normix(
    const arma::uword & k_p_, 
    const arma::uword & k_n_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_model(k_p_, k_n_, inf_, sup_, l_), 
  mean_p(arma::zeros(k_p_)), sd_p(arma::zeros(k_p_)),
  mean_n(arma::zeros(k_n_)), sd_n(arma::zeros(k_n_)){};
  
  f_normix(
    const arma::vec & w_p_, 
    const arma::vec & mean_p_, 
    const arma::vec & sd_p_,
    const arma::vec & w_n_, 
    const arma::vec & mean_n_, 
    const arma::vec & sd_n_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_model(w_p_, w_n_, inf_, sup_, l_),
  mean_p(mean_p_), sd_p(sd_p_),
  mean_n(mean_n_), sd_n(sd_n_){};
  
  virtual ~f_normix(){};
  
  // --- Methods for creating benchmark
  void render_pair_from_pc(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & p, 
      const bool & is_star, 
      const bool & has_common,
      const int & maxit_0, 
      const double & eps_0);
  
  void render_pair_from_nc(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & p, 
      const bool & is_star, 
      const bool & has_common,
      const int & maxit_0, 
      const double & eps_0);
  
  void set_rand_par_p(const arma::uword & k) {
    sd_p.at(k) = sd_n.max() * (1. + 0.5 * arma::randu());
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
  
  void set_alt_param(){};
  
  // --- Sampling methods
  double rand_1_from_pc(const arma::uword & z) {
    return mean_p.at(z) + sd_p.at(z) * arma::randn();
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
    f_pair_norm f(
        w_pair_p.at(z)/c_norm.at(z), mean_pair_p.at(z), sd_pair_p.at(z),
        w_pair_n.at(z)/c_norm.at(z), mean_pair_n.at(z), sd_pair_n.at(z),
        0., 1., lambda, 1
    );
    return sample_from_pair(f, delta, eps_d, p_accept, is_build,
                            grid, h, p_cell, counter,
                            optim, use_mono, n_points, 
                            maxit, eps_f, eps_g);
  };
  
  // --- Negative mixture CDF
  arma::vec cdf(const arma::vec & x) {
    return pnormix(x, w_p, mean_p, sd_p) - pnormix(x, w_n, mean_n, sd_n);
  }
  
  // --- Negative mixture quantile function
  arma::vec inv_cdf_init(const unsigned & n_bins) {
    return arma::linspace(arma::min(mean_p - 2.0 * sd_p), 
                          arma::max(mean_p + 2.0 * sd_p), 
                          n_bins);
  };
  
  // --- Negative mixture PDF
  double w_pdf_p(const double & x) {
    return arma::sum(w_p % arma::normpdf(x, mean_p, sd_p));
  };
  
  double w_pdf_n(const double & x) {
    return arma::sum(w_n % arma::normpdf(x, mean_n, sd_n));
  };
  
  // --- Methods for pairs
  bool is_valid_pair(const arma::uword & i, const arma::uword & j) {
    return (sd_n.at(j) < sd_p.at(i));
  };
  
  double a_star_pair(const arma::uword & i, const arma::uword & j) {
    return a_star_norm(mean_p.at(i), sd_p.at(i), mean_n.at(j), sd_n.at(j), false);
  };
  
  void set_pairs_param() {
    const arma::uvec r_0 = list_pair.row(0).t();
    const arma::uvec r_1 = list_pair.row(1).t();
    // --- Extracting mean and sd
    mean_pair_p = mean_p.elem(r_0);
    sd_pair_p = sd_p.elem(r_0);
    mean_pair_n = mean_n.elem(r_1);
    sd_pair_n = sd_n.elem(r_1);
    
    // --- Cleaning positive and negative remainder
    mean_r_p = mean_p.elem(list_r_p);
    sd_r_p = sd_p.elem(list_r_p);
    
    mean_r_n = mean_n.elem(list_r_n);
    sd_r_n = sd_n.elem(list_r_n);
  };
  
  // --- Methods for stratified sampling
  double stratified_accept_ratio(const double & y) {
    double t_prop = this->w_pdf_p(y) - sum(w_pair_n % arma::normpdf(y, mean_pair_n, sd_pair_n));
    double t_n_r = sum(r_n % arma::normpdf(y, mean_r_n, sd_r_n));
    return 1. - t_n_r / t_prop;
  };
  
  // ---
  double f_grad(Constvec & x, Refvec grad) {return 0.;};
  
public:
  arma::vec mean_p, sd_p, mean_n, sd_n;
  arma::vec mean_pair_p, sd_pair_p, mean_pair_n, sd_pair_n;
  arma::vec mean_r_p, sd_r_p;
  arma::vec mean_r_n, sd_r_n;
};

#endif /* INST_INCLUDE_F_NORMIX_HPP_ */
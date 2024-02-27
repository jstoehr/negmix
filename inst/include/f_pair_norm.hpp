#ifndef INST_INCLUDE_F_PAIR_NORM_HPP_
#define INST_INCLUDE_F_PAIR_NORM_HPP_

#include <f_pair.hpp>
#include <f_a_norm.hpp>
#include <a_star.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_pair_norm: public f_pair {
public:
  f_pair_norm(
    const double & w_p_, 
    const double & m_p, 
    const double & s_p,
    const double & w_n_, 
    const double & m_n, 
    const double & s_n,
    const double & inf_, 
    const double & sup_, 
    const double & l_, 
    const double & scale_
  ): f_pair(w_p_, w_n_, inf_, sup_, l_, scale_), 
  mean_p(m_p), sd_p(s_p), v_p(s_p * s_p),
  mean_n(m_n), sd_n(s_n), v_n(s_n * s_n){};
  
  virtual~f_pair_norm(){};
  
  // ---
  double rand_1() {
    return mean_p + sd_p * arma::randn();
  };
  
  // --- cdf function
  double w_cdf_p(const double & x) {
    return w_p * arma::normcdf(x, mean_p, sd_p);
  }
  
  double w_cdf_n(const double & x) {
    return w_n * arma::normcdf(x, mean_n, sd_n);
  }
  
  arma::vec cdf(const arma::vec & x) {
    return w_p * arma::normcdf(x, mean_p, sd_p) - w_n * arma::normcdf(x, mean_n, sd_n);
  };
  
  // --- 
  double q_p(const double & x) {
    return R::qnorm(x, mean_p, sd_p, true, false);
  }
  
  // --- pdf function
  double w_pdf_p(const double & x) {
    return w_p * arma::normpdf(x, mean_p, sd_p);
  };
  
  double w_pdf_n(const double & x) {
    return w_n * arma::normpdf(x, mean_n, sd_n);
  };
  
  // --- gradient function
  arma::vec log_pdf(const double & x) {
    return { arma::log_normpdf(x, mean_p, sd_p), arma::log_normpdf(x, mean_n, sd_n) };
  };
  
  arma::vec grad_stat(const double & x) {
    return { (mean_p - x) / v_p, (mean_n - x) / v_n };
  };
  
  double grad(const double & x) {
    return (x - mean_n) / v_n * this->w_pdf_n(x) - (x - mean_p) / v_p * this->w_pdf_p(x);
  };
  
  // --- Methods for building histogram
  void init_grid(
      const double & alpha, 
      arma::vec & base, 
      arma::ivec & mono,
      const bool & optim,
      const int & maxit, 
      const double & eps_f, 
      const double & eps_g
  );
  
  void set_step_size(
      const double & n_points, 
      arma::vec & base, 
      arma::ivec & mono,
      arma::vec & eps
  ) {
    // --- Setting the step size but also introduce points
    // --- for unbounded densities
    arma::uword n = mono.n_rows;
    double eps_ = (base.at(n) - base.at(0)) / n_points;
    eps.set_size(n);
    eps.fill(eps_);
  };
  
  void get_mono_comp(
      const arma::vec & base, 
      arma::ivec & mono_p, 
      arma::ivec & mono_n
  ) {
    set_mono_comp(base, mean_p, mono_p);
    set_mono_comp(base, mean_n, mono_n);
  };
  
  void join_grids(
      arma::field<arma::vec> & g_loc, 
      arma::field<arma::vec> & h_loc,
      arma::vec & grid, 
      arma::vec & h
  );
  
public:
  double mean_p, sd_p, v_p, mean_n, sd_n, v_n;
};

#endif /* INST_INCLUDE_F_PAIR_NORM_HPP_ */
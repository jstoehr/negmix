#ifndef INST_INCLUDE_F_PAIR_GAM_HPP_
#define INST_INCLUDE_F_PAIR_GAM_HPP_

#include <f_a_gam.hpp>
#include <a_star.hpp>
#include <dnegmix.hpp>
#include <build_grid.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_pair_gam: public f_pair {
public:
  f_pair_gam(
    const double & w_p_, 
    const double & a_p, 
    const double & b_p,
    const double & w_n_, 
    const double & a_n, 
    const double & b_n,
    const double & inf_, 
    const double & sup_, 
    const double & l_, 
    const double & scale_): f_pair(w_p_, w_n_, inf_, sup_, l_, scale_),
    alpha_p(a_p), beta_p(b_p), alpha_n(a_n), beta_n(b_n), inv_beta_p(1./b_p), inv_beta_n(1./b_n)
  {
    m_p = (alpha_p - 1)/beta_p;
    if (beta_n * (1 - alpha_p) == beta_p * (1 - alpha_n)) {
      // --- to avoid rouding issues
      m_n = m_p;
    } else {
      m_n = std::max(0., (alpha_n - 1)/beta_n);
    }
  };
  
  virtual~f_pair_gam(){};
  
  // ---
  double rand_1() {
    return arma::randg(arma::distr_param(alpha_p, inv_beta_p));
  };
  
  // --- cdf function
  double w_cdf_p(const double & x) {
    return w_p * R::pgamma(x, alpha_p, inv_beta_p, true, false);
  }
  
  double w_cdf_n(const double & x) {
    return w_n * R::pgamma(x, alpha_n, inv_beta_n, true, false);
  }
  
  arma::vec cdf(const arma::vec & x) {
    arma::vec ans(x.n_rows);
    for (arma::uword k = 0; k < x.n_rows; k++) {
      ans.at(k) = this->w_cdf_p(x.at(k)) - this->w_cdf_n(x.at(k));
      // ans.at(k) = this->cdf_1(x.at(k));
    }
    return ans;
  };
  
  // --- 
  double q_p(const double & x) {
    return R::qgamma(x, alpha_p, inv_beta_p, true, false);
  }
  
  // --- pdf function
  double w_pdf_p(const double & x) {
    return w_p * exp(log_gampdf_1(x, alpha_p, beta_p));
  };
  
  double w_pdf_n(const double & x) {
    return w_n * exp(log_gampdf_1(x, alpha_n, beta_n));
  };
  
  // --- grad function
  arma::vec log_pdf(const double & x) {
    return { log_gampdf_1(x, alpha_p, beta_p), log_gampdf_1(x, alpha_n, beta_n) };
  };
  
  arma::vec grad_stat(const double & x) {
    return { (alpha_p - 1.)/x - beta_p, (alpha_n - 1.)/x - beta_n };
  };
  
  double grad(const double & x) {
    return (this->w_pdf_p(x) * (alpha_p - 1. - beta_p * x) -
            this->w_pdf_n(x) * (alpha_n - 1. - beta_n * x)) / x;
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
  
  void insert_points(
      arma::vec & base,
      arma::ivec & mono,
      arma::vec & eps
  );
  
  void set_step_size(
      const double & n_points, 
      arma::vec & base, 
      arma::ivec & mono, 
      arma::vec & eps
  ) {
    arma::uword n = mono.n_rows;
    double eps_ = (base.at(n) - base.at(0)) / n_points;
    eps.set_size(n);
    eps.fill(eps_);
    if (alpha_p < 1.) {
      this->insert_points(base, mono, eps);
    }
  };
  
  void get_mono_comp(
      const arma::vec & base, 
      arma::ivec & mono_p, 
      arma::ivec & mono_n
  ) {
    set_mono_comp(base, m_p, mono_p);
    set_mono_comp(base, m_n, mono_n);
  };
  
  void join_grids(
      arma::field<arma::vec> & g_loc, 
      arma::field<arma::vec> & h_loc,
      arma::vec & grid, 
      arma::vec & h
  ) {
    merge_grids(g_loc, h_loc, grid, h);
  };
  
public:
  double alpha_p, beta_p, alpha_n, beta_n, inv_beta_p, inv_beta_n, m_p, m_n;
};

#endif /* INST_INCLUDE_F_PAIR_GAM_HPP_ */
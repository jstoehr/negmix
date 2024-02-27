// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_pair_gam.hpp>
#include <build_grid.hpp>
#include <optim.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void f_pair_gam::init_grid(const double & alpha, 
                           arma::vec & base, 
                           arma::ivec & mono,
                           const bool & optim, 
                           const int & maxit, 
                           const double & eps_f, 
                           const double & eps_g) {
  // --- Output on a
  double a_star = a_star_gam(alpha_p, beta_p, alpha_n, beta_n, false);
  
  // --- Output on extreme points on the grid
  double x_opt = (alpha_p - alpha_n)/(beta_p - beta_n);
  double x_1, x_n;
  if (alpha_p < 1 and (x_opt != 0. or a != a_star)) {
    // --- We move away from 0
    x_1 = R::qgamma(.5 * alpha, alpha_n, inv_beta_n, true, false);
    x_n = R::qgamma(1. - .5 * alpha, alpha_n, inv_beta_n, true, false);
  } else {
    x_1 = 0.;
    x_n = R::qgamma(1. - alpha, alpha_n, inv_beta_n, true, false);
  }
  
  // --- 
  arma::vec x;
  arma::ivec config;
  if (a - a_star < 1e-14) {
    // --- Global minimum at x_opt
    get_config_modes(m_p, m_n, x_opt, true, true, x, config);
  } else if (!optim) {
    // --- Limit value for a when component shares the same mode
    double a_lim = log_gampdf_1(m_p, alpha_n, beta_n) - log_gampdf_1(m_p, alpha_p, beta_p);
    a_lim = beta_n / beta_p * exp(a_lim);
    bool cond_a = (a < a_lim);
    // --- Configuration with modes
    get_config_modes(m_p, m_n, x_opt, cond_a, false, x, config);
  } else {
    // --- Otherwise compute the limit value for a
    // --- x-axis of the tangent point gives info
    // --- on the location of the min
    
    if (alpha_p == alpha_n and alpha_p == 1.) {
      // --- Global maximum is known
      double x_m = 0.;
      double cond = a * beta_p * beta_p/(beta_n * beta_n);
      if (cond < 1.) {
        // --- Optimal initial grid for symmetric mixture
        x_m = log(cond)/(beta_p - beta_n);
      }
      x = { x_m };
      config = { 1, -1 };
    } else {
      double x_a = 0., a_lim = 0.;
      if (alpha_p == alpha_n and alpha_p < 1.) {
        // --- Polynomial giving sign of gradient of psi
        double a_ = beta_p * beta_n * (beta_p - beta_n);
        double b_ = (alpha_p - 1.) * (beta_n * beta_n - beta_p * beta_p);
        double c_ = alpha_p * (1. - alpha_p) * (beta_n - beta_p);
        // --- Positive root
        x_a = -0.5 * (b_ + sqrt(b_ * b_ - 4. * a_ * c_)) / a_; 
        a_lim = (beta_n * x_a + 1. - alpha_n) / (beta_p * x_a + 1. - alpha_p);
        a_lim *= this->w_pdf_n(x_a)/this->w_pdf_p(x_a) * a; 
      } else if (alpha_p == alpha_n and alpha_p > 1.) {
        // --- Always a gloabl maximum in this situation
        x_a = 0.;
        a_lim = a + 1.;
      } else if (m_p == m_n) {
        x_a = m_p;
        a_lim = log_gampdf_1(m_p, alpha_n, beta_n) - log_gampdf_1(m_p, alpha_p, beta_p);
        a_lim = beta_n / beta_p * exp(a_lim);
      } else if (m_p != m_n) {
        // --- Polynomial giving sign of gradient of psi
        double a_ = beta_n * (beta_p - beta_n);
        double b_ = beta_p * (1. - alpha_n) + beta_n * (2. * alpha_n - alpha_p);
        double c_ = (alpha_n - alpha_p) * (1. - alpha_n);
        double d = sqrt(b_ * b_ - 4. * a_ * c_); 
        double inf_a, sup_a;
        // --- affine transformation is 0 at m_p
        // --- psi is negative before m_n
        if (m_p < m_n) {
          // --- sup of obj function is in (m_n, argmax psi)
          inf_a = m_n;
          sup_a = -0.5 * (b_ + d)/a_;
        } else {
          // --- sup of obj function is in (argmin psi, m_n)
          inf_a = 0.5 * (d - b_)/a_;
          sup_a = m_n;
        }
        
        f_a_gam f_a(alpha_p, beta_p, alpha_n, beta_n, inf_a, sup_a, lambda);
        x_a = 0.5 * (inf_a + sup_a);
        int res = bfgs_box(f_a, x_a, a_lim, maxit, eps_f, eps_g);
        // --- WARNING bfgs is used with minimization version
        a_lim = -a_lim;
      }
      bool cond_a = a < a_lim;
      get_config_optim(m_p, m_n, x_opt, x_a, cond_a, x, config);
    }
  }
  set_init_grid(x_1, x_n, x, config, base, mono);
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_pair_gam::insert_points(arma::vec & base, arma::ivec & mono, arma::vec & eps) {
  // --- Introducing intermediate points to avoid numeric issues
  double m_temp = this->operator()(base.at(0));
  double eps_temp = .01 / m_temp;
  double x_temp = base.at(0) + 100 * eps_temp;
  arma::uword n = 0;
  arma::vec v_eps_;
  
  if (x_temp < base.at(1)) {
    while (x_temp < base.at(1)) {
      ++n;
      v_eps_.resize(n);
      v_eps_.tail(1) = eps_temp;
      m_temp = this->operator()(x_temp);
      eps_temp = .01 / m_temp;
      x_temp += 100 * eps_temp;
    }
    
    arma::vec base_ = base.at(0) + 100 * v_eps_;
    arma::ivec mono_(n);
    mono_.fill(mono.at(0));
    
    base.insert_rows(1, base_);
    mono.insert_rows(0, mono_);
    eps.insert_rows(0, v_eps_);
  }
}

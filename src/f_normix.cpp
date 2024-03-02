// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_normix.hpp>
#include <optim.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void f_normix::render_pair_from_pc(const arma::uword & i, 
                                   const arma::uword & j, 
                                   const double & p, 
                                   const bool & is_star, 
                                   const bool & has_common,
                                   const int & maxit_0, 
                                   const double & eps_0) {
  double p_ = p;
  if (!is_star) {
    p_ = arma::randu() * p;
  }
  
  if (has_common) {
    mean_n.at(j) = mean_p.at(i);
    sd_n.at(j) = sd_p.at(i) * (1. - p_);
  } else {
    double b = arma::randu(arma::distr_param(1. - p_, 1.));
    sd_n.at(j) = sd_p.at(i) * (1. - p_)/b;
    double temp = sqrt(-log(b) * 2 * (sd_p.at(i) * sd_p.at(i) - sd_n.at(j) * sd_n.at(j)));
    if (arma::randu() < 0.5) {
      mean_n.at(j) = mean_p.at(i) - temp;
    } else {
      mean_n.at(j) = mean_p.at(i) + temp;
    }
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_normix::render_pair_from_nc(const arma::uword & i, 
                                   const arma::uword & j,
                                   const double & p, 
                                   const bool & is_star, 
                                   const bool & has_common,
                                   const int & maxit_0, 
                                   const double & eps_0) {
  double p_ = p;
  if (!is_star) {
    p_ = arma::randu() * p;
  }
  
  if (has_common) {
    mean_p.at(i) = mean_n.at(j);
    sd_p.at(i) = sd_n.at(j) / (1. - p_);
  } else {
    double b = arma::randu(arma::distr_param(1. - p_, 1.));
    sd_p.at(i) = sd_n.at(j) * b/(1. - p_);
    double temp = sqrt(-log(b) * 2 * (sd_p.at(i) * sd_p.at(i) - sd_n.at(j) * sd_n.at(j)));
    if (arma::randu() < 0.5) {
      mean_p.at(i) = mean_n.at(j) - temp;
    } else {
      mean_p.at(i) = mean_n.at(j) + temp;
    }
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_normix::add_single_comp(const arma::uword & k_1, const arma::uword & k_2) {
  // --- Range of mean values
  double min = mean_p.head(k_p - k_1 - k_2).min();
  double max = mean_p.head(k_p - k_1 - k_2).max();
  if (min == max) {
    min = std::min(mean_n.min(), min);
    max = 1.1 * std::max(mean_n.max(), max);
  }
  
  // --- Containers
  arma::vec m_t, s_t;
  
  if (k_1 > 0) {
    m_t = arma::randu(k_1, arma::distr_param(min, max));
    arma::uvec indices = arma::randi<arma::uvec>(k_1, arma::distr_param(0, k_p - k_1 - k_2 - 1));
    s_t = sd_p.elem(indices) % arma::randu(k_1);
  }
  
  if (k_2 > 0) {
    arma::vec m_t_2 = arma::randu(k_2, arma::distr_param(min, max));
    arma::vec s_t_2 = sd_n.min() * arma::randu(k_2);
    m_t = arma::join_vert(m_t, m_t_2);
    s_t = arma::join_vert(s_t, s_t_2);
  }
  mean_p.tail(k_1 + k_2) = m_t;
  sd_p.tail(k_1 + k_2) = s_t;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


arma::vec f_normix::w_valid_mixt(const arma::uword & i, 
                                 const arma::uword & j, 
                                 const double & w_p_loc, 
                                 const double & w_n_loc,
                                 const int & maxit_0, 
                                 const double & eps_0, 
                                 const int & maxit, 
                                 const double & eps_f, 
                                 const double & eps_g) {
  double m_p = mean_p.at(i), s_p = sd_p.at(i);
  double m_n = mean_n.at(j), s_n = sd_n.at(j);
  // --- Param of all other positive weight component
  arma::vec w_p_ = 1. + arma::randu(k_p - 1);
  arma::vec mean_p_ = mean_p;
  arma::vec sd_p_ = sd_p;
  mean_p_.shed_row(i);
  sd_p_.shed_row(i);
  
  // --- Interval where the two components signed mixture is negative
  double ivp = 1. / (s_p * s_p);
  double ivn = 1. / (s_n * s_n);
  double a = 0.5 * (ivn - ivp);
  double b = m_p * ivp - m_n * ivn;
  double c = 0.5 * (ivn * m_n * m_n - ivp * m_p * m_p) + log(s_n/s_p) + log(w_p_loc/w_n_loc);
  double d = b * b - 4 * a * c;
  double inf = (-sqrt(d) - b) / (2. * a);
  double sup = (sqrt(d) - b) / (2. * a);
  
  // --- Minimum for the pair
  f_pair_norm f2(w_n_loc, m_n, s_n, w_p_loc, m_p, s_p, inf, sup, lambda, -1);
  double x_min = 0.5 * (inf + sup);
  double f_min;
  int res = bfgs_box(f2, x_min, f_min, maxit, eps_f, eps_g);
  if (res != 0) {
    x_min = 0.5 * (inf + sup);
    // Rcpp::warning(" Failed to find the minimum on the negative support. Correction applied.");
  }
  
  // --- Optimization
  f_obj_norm f({ w_p_loc, -w_n_loc }, { m_p, m_n }, { s_p, s_n }, 
               w_p_, mean_p_, sd_p_,
               inf, sup, lambda);
  
  double x_opt = x_min, f_opt;
  res = bfgs_box(f, x_opt, f_opt, maxit, eps_f, eps_g);
  if (res != 0) {
    // Rcpp::warning(" Failed to find proper weight correction. Setting new reference weight.");
  }
  
  arma::vec w_insert;
  if (isnan(exp(-f_opt)) or isinf(exp(-f_opt)) or res != 0) {
    // --- Optimisation failed or has no solution
    w_p_.zeros();
    w_insert = { a_star_norm(m_p, s_p, m_n, s_n, false) * w_n_loc };
  } else {
    // --- Update weight of the positive weight components
    // --- in order to balance the negative one
    w_p_ *= exp(-f_opt);
    w_insert = { w_p_loc };
  }
  w_p_.insert_rows(i, w_insert);
  
  return w_p_;
}

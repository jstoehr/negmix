// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <math.h> 
#include <f_gammix.hpp>
#include <f_zero.hpp>
#include <f_obj_gam.hpp>
#include <optim.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void f_gammix::render_pair_from_pc(const arma::uword & i, 
                                   const arma::uword & j,
                                   const double & p, 
                                   const bool & is_star, 
                                   const bool & has_common,
                                   const int & maxit_0, 
                                   const double & eps_0) {
  // --- Containers
  double b_n, f, df;
  bool out;
  // --- Set a target to get acceptance rate lower than p
  double p_ = arma::randu() * p;
  double target = log(1./(1. - p_));
  // --- Find alpha_n for an a_star associated to the latter acceptance rate
  if (has_common) {
    alpha_n.at(j) = alpha_p.at(i);
    b_n = beta_p.at(i) * exp(-log(1 - p_)/alpha_p.at(i));
  } else {
    render_alpha_n f_an(alpha_p.at(i), beta_p.at(i), target, 
                        std::max(lambda, alpha_p.at(i) * 1.1 / log(2)));
    alpha_n.at(j) = alpha_p.at(i) * 1.1;
    out = newton_raphson_right(f_an, alpha_n.at(j), f, df, maxit_0, eps_0);
    b_n = alpha_n.at(j)/alpha_p.at(i) * beta_p.at(i);
    if (out or isinf(alpha_n.at(j))) {
      alpha_n.at(j) = alpha_p.at(i) * (1. + 0.1 * arma::randu());
      b_n = beta_p.at(i) * exp(-log(1 - p_)/alpha_n.at(j));
    }
  }
  
  // --- Update beta_n to have an acceptance rate of p
  // --- (a_star is currently minimal, so a value exists)
  if (is_star) {
    p_ = p;
  } else {
    p_ = arma::randu(arma::distr_param(p_, p));
  }
  target = log(1./(1. - p_));
  render_beta_n f_bn(alpha_p.at(i), beta_p.at(i), alpha_n.at(j), target, 
                     std::max(lambda, b_n * 1.1 / log(2)));
  b_n *= 1.1;
  out = newton_raphson_right(f_bn, b_n, f, df, maxit_0, eps_0);
  beta_n.at(j) = b_n;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_gammix::render_pair_from_nc(const arma::uword & i, 
                                   const arma::uword & j,
                                   const double & p, 
                                   const bool & is_star, 
                                   const bool & has_common,
                                   const int & maxit_0, 
                                   const double & eps_0) {
  // --- Containers
  double b_p, f, df;
  bool out;
  // --- Set a target to get acceptance rate lower than p
  double p_ = arma::randu() * p;
  double target = log(1./(1. - p_));
  // --- Find alpha_n for an a_star associated to the latter acceptance rate
  if (has_common) {
    alpha_p.at(i) = alpha_n.at(j);
    b_p = beta_n.at(j) * exp(log(1 - p_)/alpha_p.at(i));
  } else {
    render_alpha_p f_ap(alpha_n.at(j), beta_n.at(j), target, lambda);
    out = newton_raphson_box(f_ap, alpha_p.at(i), f, df, maxit_0, eps_0);
    b_p = alpha_p.at(i)/alpha_n.at(j) * beta_n.at(j);
    if (out) {
      // --- Equality case
      b_p = beta_n.at(j) * exp(log(1 - p_)/alpha_p.at(i));
    }
  }
  
  // --- Update beta_n to have an acceptance rate of p
  // --- (a_star is currently minimal, so a value exists)
  if (is_star) {
    p_ = p;
  } else {
    p_ = arma::randu(arma::distr_param(p_, p));
  }
  target = log(1./(1. - p_));
  render_beta_p f_bp(alpha_p.at(i), alpha_n.at(j), beta_n.at(j), target, lambda);
  beta_p.at(i) = b_p;
  out = newton_raphson_box(f_bp, b_p, f, df, maxit_0, eps_0);
  if (!out) {
    beta_p.at(i) = b_p;
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_gammix::add_single_comp(const arma::uword & k_1, const arma::uword & k_2) {
  // --- Range of alpha values
  double min = alpha_p.head(k_p - k_1 -k_2).min();
  double max = alpha_p.max();
  if (min == max) {
    // --- When only one positive weight component is involved in all pairs
    min = 0.;
  }
  
  // --- Containers
  arma::vec a_t(k_1), b_t(k_1);
  
  if (k_1 > 0) {
    a_t = arma::randu(k_1, arma::distr_param(min, max));
    arma::vec prob;
    arma::uvec indices;
    for (arma::uword k = 0; k < k_1; k++) {
      indices = find(alpha_n >= a_t.at(k));
      prob.set_size(indices.n_rows);
      prob.fill(1./indices.n_rows);
      b_t.at(k) = beta_n.at(sample_int_1(indices, prob)) * arma::randu();
    }
  }
  
  if (k_2 > 0) {
    arma::vec a_t_2 = arma::randu(k_2, arma::distr_param(min, alpha_n.max()));
    arma::vec b_t_2(k_2);
    arma::uvec indices;
    for (arma::uword k = 0; k < k_2; k++) {
      indices = find(alpha_n >= a_t_2.at(k));
      b_t_2.at(k) = beta_n.elem(indices).max() * (1. + 0.5 * arma::randu());
    }
    a_t = join_vert(a_t, a_t_2);
    b_t = join_vert(b_t, b_t_2);
  }
  alpha_p.tail(k_1 + k_2) = a_t;
  beta_p.tail(k_1 + k_2) = b_t;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


arma::vec f_gammix::w_valid_mixt(const arma::uword & i, 
                                 const arma::uword & j, 
                                 const double & w_p_loc, 
                                 const double & w_n_loc,
                                 const int & maxit_0, 
                                 const double & eps_0, 
                                 const int & maxit, 
                                 const double & eps_f, 
                                 const double & eps_g) {
  double a_p = alpha_p.at(i), b_p = beta_p.at(i);
  double a_n = alpha_n.at(j), b_n = beta_n.at(j);
  // --- Param of all other positive weight component
  arma::vec w_p_ = 1. + arma::randu(k_p - 1);
  arma::vec alpha_p_ = alpha_p;
  arma::vec beta_p_ = beta_p;
  alpha_p_.shed_row(i);
  beta_p_.shed_row(i);
  
  // --- Interval where the two components signed mixture is negative
  double x_star = x_star_gam(a_p, b_p, a_n, b_n);
  double inf, sup;
  
  // printf(" Search for the support of the negative part |");
  if (x_star > 0.) {
    double g_val, dg_val;
    neg_part_gam g(w_p_loc, a_p, b_p, w_n_loc, a_n, b_n, lambda);
    // --- First zero
    // printf(" First zero");
    bool out = newton_raphson_box(g, inf, g_val, dg_val, maxit_0, eps_0);
    // printf(" found | ");
    // --- Second zero
    g.left = false;
    g.inf = x_star;
    double eps = 0.1;
    double d = lambda * eps * (lambda * eps + 4. * a_n - 4. * a_p);
    sup = x_star + 0.5 * (lambda * eps + sqrt(d))/(b_n - b_p);
    // printf(" Second zero");
    out = newton_raphson_right(g, sup, g_val, dg_val, maxit_0, eps_0);
    // printf(" found |");
  } else {
    inf = 0.;
    sup = (a_n * log(b_n) - a_p * log(b_p) + log(w_n_loc) - log(w_p_loc) - eps_0)/(b_n - b_p);
  }
  
  // --- Minimum for the pair
  f_pair_gam f2(w_n_loc, a_n, b_n, w_p_loc, a_p, b_p, inf, sup, lambda, -1);
  double x_min = 0.5 * (inf + sup);
  double f_min;
  int res = bfgs_box(f2, x_min, f_min, maxit, eps_f, eps_g);
  if (res != 0) {
    x_min = 0.5 * (inf + sup);
    // Rcpp::warning(" Failed to find the minimum on the negative support. Correction applied.");
  }

  // --- Optimization
  f_obj_gam f({ w_p_loc, -w_n_loc }, { a_p, a_n }, { b_p, b_n }, 
              w_p_, alpha_p_, beta_p_,
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
    w_insert = { a_star_gam(a_p, b_p, a_n, b_n, false) * w_n_loc };
  } else {
    // --- Update weight of the positive weight components
    // --- in order to balance the negative one
    w_p_ *= exp(-f_opt);
    w_insert = { w_p_loc };
  }
  w_p_.insert_rows(i, w_insert);

  return w_p_;
}

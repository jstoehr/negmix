#ifndef INST_INCLUDE_F_OBJ_GAM_HPP_
#define INST_INCLUDE_F_OBJ_GAM_HPP_

#include <f_model.hpp>
#include <dnegmix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline void grad_log_dgammix(const double & x, 
                             const arma::vec & w, 
                             const arma::vec & alpha, 
                             const arma::vec & beta,
                             double & log_f, double & d_log_f) {
  arma::vec y(alpha.n_rows);
  for (arma::uword k = 0; k < alpha.n_rows; k++) {
    y.at(k) = log_gampdf_1(x, alpha.at(k), beta.at(k));
  }
  arma::vec temp;
  log_f = log_sum_exp(w, y, temp);
  d_log_f = arma::sum(((alpha - 1.)/x - beta) % temp)/arma::sum(temp);
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class f_obj_gam: public f_negmix {
public:
  f_obj_gam(
    const arma::vec & w_pair_, 
    const arma::vec & alpha_pair_, 
    const arma::vec & beta_pair_,
    const arma::vec & w_p_, 
    const arma::vec & alpha_p_, 
    const arma::vec & beta_p_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_negmix(inf_, sup_, l_),
  w_pair(w_pair_), alpha_pair(alpha_pair_), beta_pair(beta_pair_), 
  w_p(w_p_), alpha_p(alpha_p_), beta_p(beta_p_){};
  
  virtual ~f_obj_gam(){};
  
  double f_grad(Constvec& x, Refvec grad) {
    double x_ = inv_logit(x[0], inf, sup, lambda);
    if (x_ == inf) {
      x_ = inv_logit(x[0] + grad[0], inf, sup, lambda);
    }
    double log_f1, d_log_f1, log_f2, d_log_f2;
    grad_log_dgammix(x_, -w_pair, alpha_pair, beta_pair, log_f1, d_log_f1);
    grad_log_dgammix(x_, w_p, alpha_p, beta_p, log_f2, d_log_f2);
    grad[0] = d_log_f2 - d_log_f1;
    // --- Accounting for the sigmoid
    grad[0] *= grad_inv_logit(x[0], inf, sup, lambda);
    // printf("x %e y %e f1 %f f2 %f f %f df %f\n",
    //        x[0], x_, log_f1, log_f2, log_f2 - log_f1, grad[0]);
    return log_f2 - log_f1;
  }
  
public:
  arma::vec w_pair, alpha_pair, beta_pair;
  arma::vec w_p, alpha_p, beta_p;
};

#endif /* INST_INCLUDE_F_OBJ_GAM_HPP_ */
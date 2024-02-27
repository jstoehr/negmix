#ifndef INST_INCLUDE_F_OBJ_NORM_HPP_
#define INST_INCLUDE_F_OBJ_NORM_HPP_

#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline void grad_log_dnormix(const double & x, 
                             const arma::vec & w, 
                             const arma::vec & mean, 
                             const arma::vec & sd,
                             double & log_f, double & d_log_f) {
  arma::vec x_c = x - mean;
  arma::vec s2 = sd % sd;
  arma::vec y = -0.5 * (x_c % x_c) / s2;
  arma::vec temp;
  log_f = log_sum_exp(w/sd, -0.5 * (x_c % x_c) / s2, temp);
  d_log_f = arma::sum(-x_c/s2 % temp)/arma::sum(temp);
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class f_obj_norm: public f_negmix {
public:
  f_obj_norm(
    const arma::vec & w_pair_, 
    const arma::vec & mean_pair_, 
    const arma::vec & sd_pair_,
    const arma::vec & w_p_, 
    const arma::vec & mean_p_, 
    const arma::vec & sd_p_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_negmix(inf_, sup_, l_),
  w_pair(w_pair_), mean_pair(mean_pair_), sd_pair(sd_pair_), 
  w_p(w_p_), mean_p(mean_p_), sd_p(sd_p_){};
  
  virtual ~f_obj_norm(){};
  
  double f_grad(Constvec& x, Refvec grad) {
    double x_ = inv_logit(x[0], inf, sup, lambda);
    double log_f1, d_log_f1, log_f2, d_log_f2;
    grad_log_dnormix(x_, -w_pair, mean_pair, sd_pair, log_f1, d_log_f1);
    grad_log_dnormix(x_, w_p, mean_p, sd_p, log_f2, d_log_f2);
    grad[0] = d_log_f2 - d_log_f1;
    // --- Accounting for the sigmoid
    grad[0] *= grad_inv_logit(x[0], inf, sup, lambda);
    // printf("inf %f sup %f x %f y %f f1 %f f2 %f f %f df %f\n", 
    //        inf, sup, x[0], x_, log_f1, log_f2, log_f2 - log_f1, grad[0]);
    return log_f2 - log_f1;
  }
  
public:
  arma::vec w_pair, mean_pair, sd_pair;
  arma::vec w_p, mean_p, sd_p;
};

#endif /* INST_INCLUDE_F_OBJ_NORM_HPP_ */
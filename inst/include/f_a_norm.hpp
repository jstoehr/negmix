#ifndef INST_INCLUDE_F_A_NORM_HPP_
#define INST_INCLUDE_F_A_NORM_HPP_

#include <f_a.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_a_norm: public f_a {
public:
  // --- Function to minimize for monotonicity condition on a
  f_a_norm(
    const double & mu_, 
    const double & sigma_, 
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_a(inf_, sup_, l_), mu(mu_), log_s(log(sigma_)), s_2(sigma_ * sigma_){};
  
  virtual~f_a_norm(){};
  
  void f_p(const double & x, double & t, double & dt, double & d2t) {
    t = -0.5 * x * x;
    dt = -1. * x;
    d2t = -1.;
  };
  
  void f_n(const double & x, double & s, double & ds, double & d2s) {
    double temp = (x - mu);
    s = -0.5 * temp * temp / s_2 - log_s;
    ds = -1. * temp / s_2;
    d2s = -1. / s_2;
  };
  
public:
  double mu, log_s, s_2;
};

#endif /* INST_INCLUDE_F_A_NORM_HPP_ */
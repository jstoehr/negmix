#ifndef INST_INCLUDE_F_A_GAM_HPP_
#define INST_INCLUDE_F_A_GAM_HPP_

#include <f_a.hpp>
#include <dnegmix.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline void f_stat_gamma(const double & x, 
                         const double & alpha, 
                         const double & beta, 
                         double & t, 
                         double & dt, 
                         double & d2t) {
  t = log_gampdf_1(x, alpha, beta);
  double temp = (alpha - 1.) / x;
  dt = temp - beta;
  d2t = -1. * temp / x;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_a_gam: public f_a {
public:
  // --- Function to minimize for monotonicity condition on a
  
  f_a_gam(
    const double & a_p, 
    const double & b_p, 
    const double & a_n, 
    const double & b_n,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_a(inf_, sup_, l_), alpha_p(a_p), beta_p(b_p), alpha_n(a_n), beta_n(b_n){};
  
  virtual~f_a_gam(){};
  
  void f_p(const double & x, double & t, double & dt, double & d2t) {
    f_stat_gamma(x, alpha_p, beta_p, t, dt, d2t);
  };
  
  void f_n(const double & x, double & s, double & ds, double & d2s) {
    f_stat_gamma(x, alpha_n, beta_n, s, ds, d2s);
  };
  
public:
  double alpha_p, beta_p, alpha_n, beta_n;
};

#endif /* INST_INCLUDE_F_A_GAM_HPP_ */
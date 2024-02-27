#ifndef INST_INCLUDE_F_A_HPP_
#define INST_INCLUDE_F_A_HPP_

#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_a: public f_negmix {
public:
  f_a(
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_negmix(inf_, sup_, l_){};
  
  virtual ~f_a(){};
  
  // --- Function to optimize to find conditions on a
  virtual void f_p(const double & x, double & t, double & dt, double & d2t){};
  virtual void f_n(const double & x, double & s, double & ds, double & d2s){};
  
  double f_grad(Constvec& x, Refvec grad) {
    double x_ = inv_logit(x[0], inf, sup, lambda);
    double t, dt, d2t, s, ds, d2s;
    this->f_p(x_, t, dt, d2t);
    this->f_n(x_, s, ds, d2s);
    double ratio = exp(s - t);
    double t_1 = d2s * dt - ds * d2t;
    double t_2 = ds * dt * (ds - dt);
    grad[0] = -1. * (t_1 + t_2)/(dt * dt) * ratio;
    // --- Accounting for the sigmoid
    grad[0] *= grad_inv_logit(x[0], inf, sup, lambda);
    return -1. * ratio * ds / dt;
  }
};

#endif /* INST_INCLUDE_F_A_HPP_ */
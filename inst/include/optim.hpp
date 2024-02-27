#ifndef INST_INCLUDE_OPTIM_HPP_
#define INST_INCLUDE_OPTIM_HPP_

#include <f_model.hpp>
#include <f_zero.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

bool newton_raphson(f_zero & fun, double & x, double & f, double & df,
                    const int & maxit_0, const double & eps_0);

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


inline bool newton_raphson_box(f_zero & f, 
                               double & x_0,
                               double & f_0,
                               double & df_0,
                               const int & maxit_0,
                               const double & eps_0) {
  x_0 = logit(0.5 * (f.inf + f.sup), f.inf, f.sup, f.lambda);
  bool out = newton_raphson(f, x_0, f_0, df_0, maxit_0, eps_0);
  x_0 = inv_logit(x_0, f.inf, f.sup, f.lambda);
  return out;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


inline bool newton_raphson_right(f_zero & f, 
                                 double & x_0,
                                 double & f_0,
                                 double & df_0,
                                 const int & maxit_0,
                                 const double & eps_0) {
  x_0 = f.lambda * log(x_0 - f.inf);
  bool out = newton_raphson(f, x_0, f_0, df_0, maxit_0, eps_0);
  x_0 = f.inf + exp(x_0/f.lambda);
  return out;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


inline int bfgs_box(f_negmix & f, 
                    double & x_opt,
                    double & f_opt,
                    const int & maxit, 
                    const double & eps_f, 
                    const double & eps_g) {
  Eigen::VectorXd x(1);
  // x[0] = logit(0.5 * (f.inf + f.sup), f.inf, f.sup, f.lambda);
  x[0] = logit(x_opt, f.inf, f.sup, f.lambda);
  int res = Numer::optim_lbfgs(f, x, f_opt, maxit, eps_f, eps_g);
  x_opt = inv_logit(x[0], f.inf, f.sup, f.lambda);
  return res;
}

#endif /* INST_INCLUDE_OPTIM_HPP_ */
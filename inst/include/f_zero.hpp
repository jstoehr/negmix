#ifndef INST_INCLUDE_F_ZERO_HPP_
#define INST_INCLUDE_F_ZERO_HPP_

#include <dnegmix.hpp>
#include <a_star.hpp>
#include <utils.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_zero {
public:
  f_zero(
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):inf(inf_), sup(sup_), lambda(l_){};
  
  virtual ~f_zero(){};
  
  double operator()(const double & x){return 0.;};
  virtual bool grad(double & x, double & f, double & df){return false;};
  
public:
  double inf, sup, lambda;
};


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class render_alpha_n: public f_zero {
public:
  render_alpha_n(
    const double & a_p, 
    const double & b_p, 
    const double & t_,
    const double & l_
  ):f_zero(a_p, 999999, l_), alpha_p(a_p), beta_p(b_p), target(t_){};
  
  virtual ~render_alpha_n(){};
  
  bool grad(double & x, double & f, double & df) {
    // --- Constraint alpha_n in (alpha_p, +infty)
    double y = inf + exp(x/lambda);
    if (y == alpha_p) {
      return true;
    } else {
      f = a_star_gam(alpha_p, beta_p, y, y/alpha_p * beta_p, true) - target;
      df = log(y) - R::digamma(y);
      df *= grad_transform_right(x, inf, lambda);
      // printf("%f ran x %f y %f f %f df %f and %f\n", lambda, x, y, f, df, f/df);
      return false;
    }
  };
  
public:
  double alpha_p, beta_p, target;
};


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class render_alpha_p: public f_zero {
public:
  render_alpha_p(
    const double & a_n, 
    const double & b_n, 
    const double & t_, 
    const double & l_
  ):f_zero(0., a_n, l_), alpha_n(a_n), beta_n(b_n), target(t_){};
  
  virtual ~render_alpha_p(){};
  
  bool grad(double & x, double & f, double & df) {
    // --- Constraint alpha_p in (0, alpha_n)
    double y = inv_logit(x, inf, sup, lambda);
    if (y == alpha_n) {
      return true;
    } else {
      f = a_star_gam(y, y/alpha_n * beta_n, alpha_n, beta_n, true) - target;
      df = R::digamma(y) - log(y);
      df *= grad_inv_logit(x, inf, sup, lambda);
      // printf("%f rap x %f y %f f %f df %f and %f\n", lambda, x, y, f, df, f/df);
      return false;
    }
  };
  
public:
  double alpha_n, beta_n, target;
};


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class render_beta_n: public f_zero {
public:
  render_beta_n(
    const double & a_p, 
    const double & b_p, 
    const double & a_n, 
    const double & t_, 
    const double & l_
  ):f_zero(b_p, 999999, l_), alpha_p(a_p), beta_p(b_p), alpha_n(a_n), target(t_){};
  
  virtual ~render_beta_n(){};
  
  bool grad(double & x, double & f, double & df) {
    // --- Constraint beta_n in (beta_p, +infty)
    double y = inf + exp(x/lambda);
    if (y == beta_p) {
      x = x + f/df;
      return true;
    } else {
      f = a_star_gam(alpha_p, beta_p, alpha_n, y, true) - target;
      df = (alpha_p - alpha_n)/(y - beta_p) + alpha_n/y;
      df *= grad_transform_right(x, inf, lambda);
      // printf("%f rbn x %f y %f f %f df %f and %f\n", lambda, x, y, f, df, f/df);
      return false;
    }
  };
  
public:
  double alpha_p, beta_p, alpha_n, target;
};


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class render_beta_p: public f_zero {
public:
  render_beta_p(
    const double & a_p, 
    const double & a_n, 
    const double & b_n, 
    const double & t_, 
    const double & l_
  ):f_zero(0., b_n, l_), alpha_p(a_p), alpha_n(a_n), beta_n(b_n), target(t_){};
  
  virtual ~render_beta_p(){};
  
  bool grad(double & x, double & f, double & df) {
    // --- Constraint beta_p in (0, beta_n)
    double y = inv_logit(x, inf, sup, lambda);
    if (y == beta_n or y == 0.) {
      x = x + f/df;
      return true;
    } else {
      f = a_star_gam(alpha_p, y, alpha_n, beta_n, true) - target;
      df = (alpha_n - alpha_p)/(beta_n - y) - alpha_p/y;
      df *= grad_inv_logit(x, inf, sup, lambda);
      // printf("%f rbp x %f y %f f %f df %f and %f\n", lambda, x, y, f, df, f/df);
      return false;
    }
  };
  
public:
  double alpha_p, alpha_n, beta_n, target;
};


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class neg_part_gam: public f_zero {
public:
  neg_part_gam(
    const double & w_p_, 
    const double & a_p, 
    const double & b_p, 
    const double & w_n_, 
    const double & a_n, 
    const double & b_n,
    const double & l_
  ):f_zero(0., (a_n - a_p)/(b_n - b_p), l_), w_p(w_p_), alpha_p(a_p), beta_p(b_p), 
  w_n(w_n_), alpha_n(a_n), beta_n(b_n), left(true){};
  
  virtual ~neg_part_gam(){};
  
  bool grad(double & x, double & f, double & df) {
    if (left) {
      double y = inv_logit(x, inf, sup, lambda);
      if (y == 0.) {
        // --- Take a point between 0. and the previous y
        y = 0.5 * inv_logit(x + f/df, inf, sup, lambda);
        // --- Assocciated value of x
        x = logit(y, inf, sup, lambda);
        f = log_gampdf_1(y, alpha_p, beta_p) - log_gampdf_1(y, alpha_n, beta_n) + log(w_p) - log(w_n);
        df = beta_n - beta_p + (alpha_p - alpha_n)/y;
        df *= grad_inv_logit(x, inf, sup, lambda);
        // printf("[%f, %f], Issue left x %f y %f f %f df %f and %f\n", inf, sup, x, y, f, df, f/df);
        return false;
      } else {
        f = log_gampdf_1(y, alpha_p, beta_p) - log_gampdf_1(y, alpha_n, beta_n) + log(w_p) - log(w_n);
        df = beta_n - beta_p + (alpha_p - alpha_n)/y;
        df *= grad_inv_logit(x, inf, sup, lambda);
        // printf("[%f, %f], left x %f y %f f %f df %f and %f\n", inf, sup, x, y, f, df, f/df);
        return false;
      }
    } else {
      double y = inf + exp(x/lambda);
      f = log_gampdf_1(y, alpha_p, beta_p) - log_gampdf_1(y, alpha_n, beta_n) + log(w_p) - log(w_n);
      df = beta_n - beta_p + (alpha_p - alpha_n)/y;
      df *= grad_transform_right(x, inf, lambda);
      // printf("right x %f y %f f %f df %f and %f\n", x, y, f, df, f/df);
      return false;
    }
  };
  
public:
  double w_p, alpha_p, beta_p, w_n, alpha_n, beta_n;
  bool left;
};


#endif /* INST_INCLUDE_F_ZERO_HPP_ */
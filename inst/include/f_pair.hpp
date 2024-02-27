#ifndef INST_INCLUDE_F_PAIR_HPP_
#define INST_INCLUDE_F_PAIR_HPP_

#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_pair: public f_negmix {
public:
  f_pair(
    const double & w_p_, 
    const double & w_n_, 
    const double & inf_, 
    const double & sup_, 
    const double & l_, 
    const double & scale_
  ):f_negmix(inf_, sup_, l_), w_p(w_p_), w_n(w_n_), a(w_p_/w_n_), scale(scale_){};
  
  virtual ~f_pair(){};
  
  // ---
  virtual double rand_1() {return 0.;};
  
  // --- cdf function for the pair
  virtual double w_cdf_p(const double & x) {return 0.;};
  virtual double w_cdf_n(const double & x) {return 0.;};
  virtual arma::vec cdf(const arma::vec & x) {return { 0. };};
  
  // --- quantile function
  virtual double q_p(const double & x) {return 0.;};
  
  // --- pdf function for the pair
  virtual double w_pdf_p(const double & x) {return 0.;};
  virtual double w_pdf_n(const double & x) {return 0.;};
  
  double operator()(const double & x){return this->w_pdf_p(x) - this->w_pdf_n(x);};
  
  // --- gradient function for optimization
  virtual arma::vec log_pdf(const double & x) {return arma::zeros(0);};
  virtual arma::vec grad_stat(const double & x) {return arma::zeros(0);};
  
  void grad_log(
      const double & x, 
      double & log_f, 
      double & d_log_f) {
    arma::vec y = this->log_pdf(x);
    arma::vec temp;
    log_f = log_sum_exp({ w_p, -w_n }, y, temp);
    d_log_f = arma::sum(this->grad_stat(x) % temp)/arma::sum(temp);
  }
  
  // double f_grad(Constvec & x, Refvec grad) {
  //   double x_ = inv_logit(x[0], inf, sup, lambda);
  //   double log_f;
  //   this->grad_log(x_, log_f, grad[0]);
  //   grad[0] *= scale * grad_inv_logit(x[0], inf, sup, lambda);
  //   return scale * log_f;
  // }
  
  virtual double grad(const double & x) {return 0.;};
  
  double f_grad(Constvec & x, Refvec grad) {
    double x_ = inv_logit(x[0], inf, sup, lambda);
    double f = this->operator()(x_);
    grad[0] = scale * this->grad(x_) / f;
    grad[0] *= grad_inv_logit(x[0], inf, sup, lambda);
    return scale * log(f);
  }
  
  // --- Methods for building histogram
  virtual void init_grid(
      const double & alpha, 
      arma::vec & base, 
      arma::ivec & mono,
      const bool & optim, 
      const int & maxit, 
      const double & eps_f, 
      const double & eps_g
  ){};
  
  void check_mono(
      const double & x_1, 
      const double & x_2,
      const double & eps_g, 
      arma::sword & type
  );
  
  virtual void set_step_size(
      const double & n_points, 
      arma::vec & base, 
      arma::ivec & mono, 
      arma::vec & eps
  ){};
  
  virtual void get_mono_comp(
      const arma::vec & base, 
      arma::ivec & mono_p, 
      arma::ivec & mono_n
  ){};
  
  double diff_h_1(
      const double & a, 
      const double & b,
      const arma::sword & mono_p, 
      const arma::sword & mono_n
  );
  
  void diff_h(
      const double & a, 
      const double & b, 
      const double & eps,
      const arma::sword & mono_p, 
      const arma::sword & mono_n,
      arma::vec & g, 
      arma::vec & h
  );
  
  void h_mono(
      const double & a, 
      const double & b, 
      const double & eps,
      const arma::sword & is_inc,
      arma::vec & g, 
      arma::vec & h
  );
  
  void h_single(
      const double & a, 
      const double & b, 
      const double & eps,
      const arma::sword & mono_p, 
      const arma::sword & mono_n,
      const unsigned & type,
      arma::vec & g, 
      arma::vec & h
  );
  
  virtual void join_grids(
      arma::field<arma::vec> & g_loc, 
      arma::field<arma::vec> & h_loc,
      arma::vec & grid, 
      arma::vec & h
  ){};
  
public:
  // --- scale = -1 for maximization
  double w_p, w_n, a, scale;
};


#endif /* INST_INCLUDE_F_PAIR_HPP_ */
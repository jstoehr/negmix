#ifndef INST_INCLUDE_F_MODEL_HPP_
#define INST_INCLUDE_F_MODEL_HPP_

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
#include <Eigen/Core>
#include <iostream>
#include <utils.hpp>

// Reference to a vector
typedef Eigen::Ref<Eigen::VectorXd>             Refvec;
typedef const Eigen::Ref<const Eigen::VectorXd> Constvec;

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class f_negmix: public Numer::MFuncGrad {
public:
  f_negmix(){};
  
  f_negmix(
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):inf(inf_), sup(sup_), lambda(l_){};
  
  virtual ~f_negmix(){};
  
public:
  double inf, sup, lambda;
};


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


class f_model: public f_negmix {
public:
  f_model(
    const arma::uword & k_p_, 
    const arma::uword & k_n_,
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_negmix(inf_, sup_, l_), k_p(k_p_), k_n(k_n_), w_p(arma::zeros(k_p_)), w_n(arma::zeros(k_n_)){};
  
  f_model(
    const arma::vec & w_p_, 
    const arma::vec & w_n_, 
    const double & inf_, 
    const double & sup_, 
    const double & l_
  ):f_negmix(inf_, sup_, l_), k_p(w_p_.n_rows), k_n(w_n_.n_rows), w_p(w_p_), w_n(w_n_){};
  
  virtual ~f_model(){};
  
  // --- General methods
  void normalize_weight(){
    double temp = arma::sum(w_p) - arma::sum(w_n);
    w_p /= temp;
    w_n /= temp;
  }
  
  virtual void set_alt_param(){};
  
  // --- Methods for creating benchmark
  // ------ Initializing parameter
  virtual void render_pair_from_pc(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & p, 
      const bool & is_star, 
      const bool & has_common,
      const int & maxit_0, 
      const double & eps_0
  ){};
  
  virtual void render_pair_from_nc(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & p, 
      const bool & is_star, 
      const bool & has_common,
      const int & maxit, 
      const double & eps_0
  ){};
  // ------ Adding pairs
  virtual void set_rand_par_p(const arma::uword & k){};
  
  arma::mat augment_n_pair(
      const arma::uword & n_pair, 
      const double & p_limit
  );
  
  virtual void add_single_comp(
      const arma::uword & k_1, 
      const arma::uword & k_2
  ){};
  
  // ------ Computing weight to balance part where pdf is negative
  virtual arma::vec w_valid_mixt(
      const arma::uword & i, 
      const arma::uword & j, 
      const double & w_p_loc, 
      const double & w_n_loc,
      const int & maxit_0, 
      const double & eps_0, 
      const int & maxit, 
      const double & eps_f, 
      const double & eps_g
  ) {
    return arma::zeros(0);
  };
  
  // ------ Adding a reference pair to get the targeted acceptance rate
  void ctrl_accept_rate(
      const double & a_ref, 
      const double & p_max
  );
  
  // --- Benchmark generator
  void create_benchmark_by_2(
      const arma::uword & k_init, 
      const arma::uword & k_pair,
      const arma::uword & k_single, 
      const arma::uword & n_pair,
      const double & p_min, 
      const double & p_max,
      const double & p_star, 
      const double & p_common,
      const double & p_rank, 
      const double & p_limit,
      const int & maxit_0, 
      const double & eps_0
  );
  
  void create_benchmark(
      const arma::uword & k_init,
      const double & p_min,
      const double & p_max,
      const double & p_star,
      const double & p_common,
      const double & p_rank,
      const int & maxit_0, 
      const double & eps_0, 
      const int & maxit, 
      const double & eps_f, 
      const double & eps_g
  );
  
  // --- Sampling methods
  virtual double rand_1_from_pc(const arma::uword & z) {return 0.;};
  
  virtual double rand_from_pair(
      const arma::uword & z, 
      const double & delta, 
      const double & eps_d, 
      const double & p_accept,
      arma::uword & is_build, 
      arma::vec & grid,
      arma::vec & h, 
      arma::vec & p_cell, 
      double & counter,
      const bool & optim, 
      const bool & use_mono,
      const double & n_points, 
      const int & maxit, 
      const double & eps_f, 
      const double & eps_g
  ) {
    return 0.;
  };
  
  // --- Negative mixture CDF
  virtual arma::vec cdf(const arma::vec & x) {return { 0. };};
  
  // --- Negative mixture quantile function
  virtual arma::vec inv_cdf_init(const unsigned & n_bins) {return { 0. };};
  
  arma::vec inv_cdf(
      const arma::vec & u, 
      const double & precision, 
      const unsigned & n_bins
  );
  
  // --- Negative mixture PDF
  virtual double w_pdf_p(const double & x) {return 0.;};
  virtual double w_pdf_n(const double & x) {return 0.;};
  
  // --- Methods for setting pairs
  virtual bool is_valid_pair(
      const arma::uword & i, 
      const arma::uword & j
  ) {return false;};
  
  virtual double a_star_pair(
      const arma::uword & i, 
      const arma::uword & j
  ) {return false;};
  
  arma::mat map_pairs(const double & p_limit);
  
  void build_pairs(
      const double & delta, 
      const double & tol_simplex
  );
  
  virtual void set_pairs_param(){};
  
  double p_accept_pairing(const double & delta) {
    arma::vec temp = (1. - delta) * w_pair_p - w_pair_n;
    return 1./(arma::sum(w_p) + arma::sum(temp.elem(find(temp < 0.)))/delta);
  };
  
  // --- Methods for stratified sampling
  virtual double stratified_accept_ratio(const double & y) {return 0.;}; 
  
public:
  arma::uword k_p, k_n;
  arma::vec w_p, w_n;
  arma::vec w_pair_p, w_pair_n, c_norm, r_p, r_n;
  arma::umat list_pair;
  arma::vec p_accept_pair;
  arma::uvec list_r_p, list_r_n;
};

#endif /* INST_INCLUDE_F_MODEL_HPP_ */
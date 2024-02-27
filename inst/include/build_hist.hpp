#ifndef INST_INCLUDE_BUILD_HIST_HPP_
#define INST_INCLUDE_BUILD_HIST_HPP_

#include <f_pair.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline void check_hist(const arma::vec & g, 
                       const arma::vec & h,
                       const double & target, 
                       const double & tol,
                       double & value, 
                       arma::uword & keep_going) {
  value = arma::sum(h % diff(g));
  if (value - target < tol) {
    keep_going = 0;
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void build_hist(f_pair & f,
                const double & delta,
                const double & eps_d,
                arma::vec & grid,
                arma::vec & h,
                arma::vec & p_cell,
                const bool & optim,
                const bool & use_mono,
                const double & n_points,
                const int & maxit, 
                const double & eps_f, 
                const double & eps_g);

#endif /* INST_INCLUDE_BUILD_HIST_HPP_ */
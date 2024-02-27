#ifndef INST_INCLUDE_SIMPLEX_HPP_
#define INST_INCLUDE_SIMPLEX_HPP_

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void one_stage_tab(const arma::vec & w_p,
                   const arma::vec & w_n,
                   const arma::umat & pair,
                   const arma::vec & a_star,
                   const double & delta,
                   arma::rowvec & x,
                   arma::mat & A,
                   arma::vec & b);

void simplex(const double & tol,
             arma::rowvec & x,
             arma::mat & A, 
             arma::vec & b,
             arma::uvec & basic_var,
             double & value);

#endif /* INST_INCLUDE_SIMPLEX_HPP_ */
#ifndef INST_INCLUDE_AR_STRATIFIED_HPP_
#define INST_INCLUDE_AR_STRATIFIED_HPP_

#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void ar_stratified(f_model & f,
                   const double & delta,
                   const double & eps_d,
                   arma::vec & ans,
                   double & p_accept,
                   double & cpu_time,
                   const bool & optim,
                   const bool & use_mono,
                   const double & n_points,
                   const int & maxit, 
                   const double & eps_f, 
                   const double & eps_g,
                   const double & tol_simplex);

#endif /* INST_INCLUDE_AR_STRATIFIED_HPP_ */
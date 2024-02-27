#ifndef INST_INCLUDE_AR_FOR_PAIR_HPP_
#define INST_INCLUDE_AR_FOR_PAIR_HPP_

#include <f_pair.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

double ar_for_pair_1(f_pair & f,
                     double & counter);

double sample_from_pair(f_pair & f,
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
                        const double & eps_g);

#endif /* INST_INCLUDE_AR_FOR_PAIR_HPP_ */
#ifndef INST_INCLUDE_BUILD_GRID_HPP_
#define INST_INCLUDE_BUILD_GRID_HPP_

#include <f_pair.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void get_config_modes(const double & m_p, 
                      const double & m_n,
                      const double & x_opt,
                      const bool & cond_a,
                      const bool & is_a_star,
                      arma::vec & x,
                      arma::ivec & config);

void get_config_optim(const double & m_p, 
                      const double & m_n,
                      const double & x_opt,
                      const double & x_a,
                      const bool & cond_a,
                      arma::vec & x,
                      arma::ivec & config);

void set_init_grid(const double & x_1, 
                   const double & x_n,
                   const arma::vec & x,
                   const arma::ivec & config,
                   arma::vec & base,
                   arma::ivec & mono);

void build_grid(f_pair & f,
                const double & delta,
                const double & epsilon,
                arma::vec & base,
                arma::ivec & mono,
                const bool & optim,
                const bool & use_mono,
                const int & maxit, 
                const double & eps_f, 
                const double & eps_g);

#endif /* INST_INCLUDE_BUILD_GRID_HPP_ */
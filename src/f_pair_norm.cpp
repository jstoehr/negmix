// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_pair_norm.hpp>
#include <build_grid.hpp>
#include <optim.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void f_pair_norm::init_grid(const double & alpha, 
                            arma::vec & base, 
                            arma::ivec & mono,
                            const bool & optim, 
                            const int & maxit, 
                            const double & eps_f, 
                            const double & eps_g) {
  // --- Output on extreme points on the grid
  double a_star = a_star_norm(mean_p, sd_p, mean_n, sd_n, false);
  double x_1 = R::qnorm(.5 * alpha, mean_n, sd_n, true, false);
  double x_n;
  if (mean_p - mean_n == 0.) {
    x_n = mean_n;
  } else {
    x_n = 2. * mean_n - x_1;
  }
  double x_opt = (a_star * v_p * mean_n - a * v_n * mean_p)  / (a_star * v_p - a * v_n);
  
  // --- 
  arma::vec x;
  arma::ivec config;
  
  if (a - a_star < 1e-14) {
    // --- Global minimum at x_opt
    get_config_modes(mean_p, mean_n, x_opt, true, true, x, config);
  } else if (!optim) {
    bool cond_a = (a < sd_p * v_p / (sd_n * v_n));
    // --- Configuration with modes
    get_config_modes(mean_p, mean_n, x_opt, cond_a, false, x, config);
  } else {
    // --- Otherwise compute the limit value for a
    // --- x-axis of the tangent point gives info
    // --- on the location of the min
    // --- Computation for the standard Normal signed mixture
    // --- Then transfor it back to original mixture
    double mu = std::abs(mean_n - mean_p)/sd_p;
    double sigma = sd_n/sd_p;
    double s_2 = v_n/v_p;
    double s_3 = s_2 * sigma;
    
    if (mu == 0.) {
      // --- Symmetry of the Normal distribution used in that situation
      double x_m = 0.;
      if (a < 1. / s_3) {
        // --- Optimal initial grid for symmetric mixture
        x_m = 2. * s_2 * log(a * s_3) / (1. - s_2);
      }
      x = { x_m };
      config = { 1, -1 };
    } else {
      // --- Limit value for a
      double x_a, a_lim;
      
      double t = 1. - s_2; 
      double sup_a = 0.5 * mu + (mu + sigma * sqrt(mu * mu * s_2 + 4. * t)) / (2. * t);
      f_a_norm f_a(mu, sigma, mu, sup_a, lambda);
      x_a = 0.5 * (mu + sup_a);
      int res = bfgs_box(f_a, x_a, a_lim, maxit, eps_f, eps_g);
      // --- WARNING bfgs is used with minimization version
      a_lim = -a_lim;
      // --- Transform x_a back to th original mixture
      if (mean_n > mean_p) {
        x_a = mean_p + sd_p * x_a;
      } else {
        x_a = mean_p - sd_p * x_a;
      }
      
      bool cond_a = a < a_lim;
      get_config_optim(mean_p, mean_n, x_opt, x_a, cond_a, x, config);
    }
  }
  set_init_grid(x_1, x_n, x, config, base, mono);
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_pair_norm::join_grids(arma::field<arma::vec> & g_loc, 
                             arma::field<arma::vec> & h_loc,
                             arma::vec & grid, 
                             arma::vec & h) {
  if (mean_p == mean_n) {
    arma::uword n_grid = g_loc.n_slices;
    arma::vec g_temp = g_loc.at(0);
    arma::vec h_temp = h_loc.at(0);
    if (n_grid > 1) {
      // --- Merging the grid
      g_temp.resize(g_temp.n_rows - 1);
      g_temp = join_vert(g_temp, g_loc.at(1));
      h_temp = join_vert(h_temp, h_loc.at(1));
    } 
    // --- Using symmetry to complete the grid
    arma::vec g_rev = -arma::reverse(g_temp);
    arma::vec h_rev = arma::reverse(h_temp);
    
    // --- Output
    g_temp.resize(g_temp.n_rows - 1);
    grid = join_vert(g_temp, g_rev);
    h = join_vert(h_temp, h_rev);
  } else {
    merge_grids(g_loc, h_loc, grid, h);
  }
}

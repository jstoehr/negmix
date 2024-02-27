// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_pair.hpp>
#include <optim.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void f_pair::check_mono(const double & x_1, const double & x_2,
                        const double & eps_g, arma::sword & type) {
  double grad_1 = this->grad(x_1);
  double grad_2 = this->grad(x_2);
  double cond = grad_1 * grad_2;
  arma::sword sign = 1;
  if (type < 0) {
    sign = -1;
  }
  
  if (std::abs(type) == 2) {
    if (sign * grad_2 > eps_g) {
      type = sign;
    } else if (sign * grad_1 < -eps_g) {
      type = -sign;
    }
  } else if (std::abs(type) == 3) {
    if (cond > 0. and sign * grad_1 < -eps_g) {
      type = -sign;
    } else if (cond < 0. and grad_1 > eps_g) {
      type = 2;
    } else if (cond < 0. and grad_1 < -eps_g) {
      type = -2;
    }
  }
}


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
                      arma::ivec & config) {
  // --- Convention 
  // --- mono = -2 (1 min in the interval), -1 (decreasing),
  // --- 0 (unknown), 1 (increasing), 2 (1 max in the interval)
  // --- (3) is for increasing, decreasing, increasing
  // --- (-3) is for decreasing, increasing, decreasing
  double cond_eq = m_p - m_n;
  
  if (cond_eq == 0.) {
    // --- Same mode
    x = { m_p };
    if (cond_a) {
      // --- One max before and after the mode
      config = { 2, 2 };
    } else {
      // --- Increasing and decreasing
      config = { 1, -1 };
    }
  } else if (cond_eq > 0.) {
    // --- mode_p > mode_n
    x = { x_opt, m_n, m_p };
    if (is_a_star) {
      config = { 2, 1, 1, 2 };
    } else {
      config = { 3, 1, 1, 2 };
    }
  } else {
    // --- mode_p > mode_n
    x = { m_p, m_n, x_opt };
    if (is_a_star) {
      config = { 2, -1, -1, 2 };
    } else {
      config = { 2, -1, -1, -3 };
    }
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void get_config_optim(const double & m_p, 
                      const double & m_n,
                      const double & x_opt,
                      const double & x_a,
                      const bool & cond_a,
                      arma::vec & x,
                      arma::ivec & config) {
  // --- Convention 
  // --- mono = -2 (1 min in the interval), -1 (decreasing),
  // --- 0 (unknown), 1 (increasing), 2 (1 max in the interval)
  // --- (3) is for increasing, decreasing, increasing
  // --- (-3) is for decreasing, increasing, decreasing
  
  // --- Keeping m_n has solely one prupose: the difference
  // --- of histograms. Otherwise, it could be removed
  
  double cond_eq = m_p - m_n;
  
  if (cond_eq == 0.) {
    // --- Same mode
    x = { m_p };
    if (cond_a) {
      // --- One max before and after the mode
      config = { 2, 2 };
    } else {
      // --- Increasing and decreasing
      config = { 1, -1 };
    }
  } else if (cond_eq > 0.) {
    // --- mode_p > mode_n
    if (cond_a) {
      x = { x_a, x_opt, m_n, m_p };
      config = { 2, -2, 1, 1, 2 };
    } else {
      x = { m_n, m_p };
      config = { 1, 1, 2 };
    }
  } else {
    // --- mode_p < mode_n
    if (cond_a) {
      x = { m_p, m_n, x_opt, x_a };
      config = { 2, -1, -1, -2, 2 };
    } else {
      x = { m_p, m_n };
      config = { 2, -1, -1 };
    }
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void set_init_grid(const double & x_1, 
                   const double & x_n,
                   const arma::vec & x,
                   const arma::ivec & config,
                   arma::vec & base,
                   arma::ivec & mono) {
  arma::uword n = x.n_rows;
  
  if (x_n < x.at(0)) {
    base = { x_1, x_n };
    mono = { config.at(0) };
  } else if (x_1 > x.at(n - 1)) {
    base = { x_1, x_n };
    mono = { config.at(n) };
  } else {
    base.set_size(n + 2);
    mono.set_size(n + 1);
    
    // --- Initialiation
    arma::uword i = 1;
    arma::uword last = n;
    base.at(0) = x_1;
    
    for (arma::uword k = 0; k < x.n_rows; k++) {
      if (x_1 < x.at(k) and x.at(k) < x_n) {
        base.at(i) = x.at(k);
        mono.at(i - 1) = config.at(k);
        ++i;
      } else if (x.at(k) >= x_n) {
        last = k;
        break;
      }
    }
    base.at(i) = x_n;
    mono.at(i - 1) = config.at(last);
    base.resize(i + 1);
    mono.resize(i);
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void build_grid(f_pair & f,
                const double & delta,
                const double & epsilon,
                arma::vec & base,
                arma::ivec & mono,
                const bool & optim,
                const bool & use_mono,
                const int & maxit, 
                const double & eps_f, 
                const double & eps_g) {
  
  double alpha = (1. - delta * (epsilon + 1.)) * (f.a - 1.) / delta;
  
  // --- Initialize the grid and storing monotonicity info
  
  f.init_grid(alpha, base, mono,
              optim, maxit, eps_f, eps_g);
  
  arma::uword n_grid = mono.n_rows;
  
  if (use_mono or optim) {
    // --- Checking if some situation we can simplify the configuration due to border conditions
    if (std::abs(mono.at(0)) > 1) {
      f.check_mono(base.at(0), base.at(1), eps_g, mono.at(0));
    }
    if (n_grid > 1) {
      if (std::abs(mono.at(n_grid - 1)) > 1) {
        f.check_mono(base.at(n_grid - 1), base.at(n_grid), eps_g, mono.at(n_grid - 1));
      }
    }
  }
  
  if (optim) {
    int res;
    double f_opt;
    arma::vec x_opt(1);
    arma::ivec m_i = { 1, -1 };
    arma::uword n_i = 0;
    arma::vec base_temp = base;
    arma::ivec mono_temp = mono;
    
    for (arma::uword k = 0; k < n_grid; k++) {
      if (mono_temp.at(k) == 2 or mono_temp.at(k) == -2) {
        f.scale = -1. * mono_temp.at(k);
        f.inf = base_temp.at(k);
        f.sup = base_temp.at(k + 1);
        x_opt.at(0) = 0.5 * (f.inf + f.sup);
        res = bfgs_box(f, x_opt.at(0), f_opt, maxit, eps_f, eps_g);
        base.insert_rows(k + 1 + n_i, x_opt);
        mono.shed_row(k + n_i);
        mono.insert_rows(k + n_i, static_cast<arma::sword>(mono_temp.at(k)/2) * m_i);
        ++n_i;
      }
    }
  }
}

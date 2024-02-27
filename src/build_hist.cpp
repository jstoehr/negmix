// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <build_hist.hpp>
#include <f_pair.hpp>
#include <build_grid.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

double f_pair::diff_h_1(const double & a, 
                        const double & b,
                        const arma::sword & mono_p, 
                        const arma::sword & mono_n) {
  double t_p, t_n;
  if (mono_p == 1) {
    t_p = this->w_pdf_p(b);
  } else {
    t_p = this->w_pdf_p(a);
  }
  if (mono_n == 1) {
    t_n = this->w_pdf_n(a);
  } else {
    t_n = this->w_pdf_n(b);
  }
  return t_p - t_n;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_pair::diff_h(const double & a, 
                    const double & b, 
                    const double & eps,
                    const arma::sword & mono_p, 
                    const arma::sword & mono_n,
                    arma::vec & g, 
                    arma::vec & h) {
  // --- Compute height of the histogram based on solely 
  // --- monotonicity of the positive and negative part
  // --- WARNING: assume that 0 in (a, b) and mu in (a, b)
  // --- is not possible
  
  // --- The grid
  grid(a, eps, b, g);
  // --- Setting container for heights
  h.zeros(g.n_rows - 1);
  
  if (mono_p == 1) {
    for (arma::uword k = 0; k < h.n_rows; k++) {
      h.at(k) = this->w_pdf_p(g.at(k + 1));
    }
  } else {
    for (arma::uword k = 0; k < h.n_rows; k++) {
      h.at(k) = this->w_pdf_p(g.at(k));
    }
  }
  
  if (mono_n == 1) {
    for (arma::uword k = 0; k < h.n_rows; k++) {
      h.at(k) -= this->w_pdf_n(g.at(k));
    }
  } else {
    for (arma::uword k = 0; k < h.n_rows; k++) {
      h.at(k) -= this->w_pdf_n(g.at(k + 1));
    }
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_pair::h_mono(const double & a, 
                    const double & b, 
                    const double & eps, 
                    const arma::sword & is_inc,
                    arma::vec & g, 
                    arma::vec & h) {
  // --- Compute height of histogram on an interval
  // --- where mixture is monotonous
  // --- (one component is increasing, the other decreasing)
  
  // --- The grid
  grid(a, eps, b, g);
  // --- Setting container for heights
  h.zeros(g.n_rows - 1);
  
  if (is_inc == 1) {
    for (arma::uword k = 0; k < h.n_rows; k++) {
      h.at(k) = this->operator()(g.at(k + 1));
    }
  } else {
    for (arma::uword k = 0; k < h.n_rows; k++) {
      h.at(k) = this->operator()(g.at(k));
    }
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void f_pair::h_single(const double & a, 
                      const double & b, 
                      const double & eps,
                      const arma::sword & mono_p, 
                      const arma::sword & mono_n,
                      const unsigned & type, 
                      arma::vec & g, 
                      arma::vec & h) {
  // --- Compute height of histogram on an interval where 
  // --- mixture has a single max (type = 1) or a single 
  // --- min (type = 0)
  
  // --- The grid
  grid(a, eps, b, g);
  // --- Setting container for heights
  h.zeros(g.n_rows - 1);
  // --- Boolean to find the max
  bool before = true;
  for (arma::uword k = 0; k < h.n_rows; k++) {
    if (before) {
      // --- Increasing part
      if ((2 * type - 1) * this->grad(g.at(k + 1)) > 0.) {
        h.at(k) = this->operator()(g.at(k + type));
      } else {
        before = false;
        h.at(k) = this->diff_h_1(g.at(k), g.at(k + 1), mono_p, mono_n);
      }
    } else {
      // --- Decreasing part
      h.at(k) = this->operator()(g.at(k + 1 - type));
    }
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
                const double & eps_g) {
  arma::vec base;
  arma::ivec mono, mono_p, mono_n;
  
  build_grid(f, delta, eps_d, base, mono, 
             optim, use_mono, maxit, eps_f, eps_g);
  
  // --- Setting step size for each subset
  arma::vec eps(mono.n_rows);
  f.set_step_size(n_points, base, mono, eps);
  
  // --- From now grid objects are not modified anymore
  arma::uword n_grid = mono.n_rows;
  double tol = eps_d / static_cast<double>(n_grid + 1);
  
  // --- Target probability on each subset
  arma::vec cdf = f.cdf(base);
  arma::vec p = arma::diff(cdf);
  
  // --- Containers
  arma::field<arma::vec> g_loc(n_grid);
  arma::field<arma::vec> h_loc(n_grid);
  arma::vec val(n_grid);
  arma::uvec keep_go = arma::ones<arma::uvec>(n_grid);
  
  double target = 1./delta - f.w_p + f.w_cdf_p(base.at(n_grid)) - f.w_cdf_p(base.at(0));
  double temp = target + 1.;
  
  //double count = 0.;
  if (optim) {
    // --- All monotonic information are used
    while (temp > target and arma::sum(keep_go) > 0) {
      // ++count;
      for (arma::uword k = 0; k < n_grid; k++) {
        if (keep_go.at(k) == 1) {
          
          f.h_mono(base.at(k), base.at(k + 1), eps.at(k), mono.at(k),
                   g_loc.at(k), h_loc.at(k));
          
          check_hist(g_loc.at(k), h_loc.at(k), p.at(k), tol,
                     val.at(k), keep_go.at(k));
          // eps.at(k) *= .5;
          eps.at(k) *= tol / (val.at(k) - p.at(k));
        }
      }
      temp = arma::sum(val);
    }
    // --- End optim case
  } else if (use_mono) {
    f.get_mono_comp(base, mono_p, mono_n);
    // --- Partial information is used
    unsigned type;
    while (temp > target and arma::sum(keep_go) > 0) {
      //++count;
      for (arma::uword k = 0; k < n_grid; k++) {
        if (std::abs(mono.at(k)) == 2 and keep_go.at(k) == 1) {
          type = (mono.at(k) + 2)/4;
          f.h_single(base.at(k), base.at(k + 1), eps.at(k),
                     mono_p.at(k), mono_n.at(k), type,
                     g_loc.at(k), h_loc.at(k));
          
          check_hist(g_loc.at(k), h_loc.at(k), p.at(k), tol,
                     val.at(k), keep_go.at(k));
          // eps.at(k) *= .5;
          eps.at(k) *= tol / (val.at(k) - p.at(k));
        } else if (std::abs(mono.at(k)) == 1 and keep_go.at(k) == 1) {
          f.h_mono(base.at(k), base.at(k + 1), eps.at(k),
                   mono.at(k),
                   g_loc.at(k), h_loc.at(k));
          
          check_hist(g_loc.at(k), h_loc.at(k), p.at(k), tol,
                     val.at(k), keep_go.at(k));
          
          // eps.at(k) *= .5;
          eps.at(k) *= tol / (val.at(k) - p.at(k));
        } else if (keep_go.at(k) == 1) {
          f.diff_h(base.at(k), base.at(k + 1), eps.at(k),
                   mono_p.at(k), mono_n.at(k), 
                   g_loc.at(k), h_loc.at(k));
          
          check_hist(g_loc.at(k), h_loc.at(k), p.at(k), tol,
                     val.at(k), keep_go.at(k));
          
          // eps.at(k) *= .5;
          eps.at(k) *= tol / (val.at(k) - p.at(k));
        }
      }
      temp = arma::sum(val);
    }
  } else {
    f.get_mono_comp(base, mono_p, mono_n);
    // --- Simple difference of the positive and negative histograms
    while (temp > target and arma::sum(keep_go) > 0) {
      // ++count;
      for (arma::uword k = 0; k < n_grid; k++) {
        if (keep_go.at(k) == 1) {
          
          f.diff_h(base.at(k), base.at(k + 1), eps.at(k),
                   mono_p.at(k), mono_n.at(k), 
                   g_loc.at(k), h_loc.at(k));
          
          check_hist(g_loc.at(k), h_loc.at(k), p.at(k), tol,
                     val.at(k), keep_go.at(k));
          // eps.at(k) *= .5;
          eps.at(k) *= tol / (val.at(k) - p.at(k));
        }
      }
      temp = arma::sum(val);
    }
  }
  
  // --- Merging grids (can be spezialised 
  // --- for instances where simplifaction 
  // --- where used, e.g., Normal signed
  // --- mixture for mean_p = mean_n)
  f.join_grids(g_loc, h_loc, grid, h);
  
  arma::vec p_temp = f.cdf(grid);
  
  p_cell.zeros(h.n_rows + 2);
  p_cell.at(0) = p_temp.at(0);
  p_cell.subvec(1, h.n_rows) = arma::diff(p_temp);
  p_cell.at(p_cell.n_rows - 1) = 1. - p_temp.at(p_temp.n_rows - 1);
}

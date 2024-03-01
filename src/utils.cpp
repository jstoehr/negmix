// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

arma::uword sample_int_1(const arma::uvec & index,
                         const arma::vec & prob) {
  double u = arma::randu();
  double p = prob.at(0);
  arma::uword i = 0;
  while (u > p) {
    i++;
    p += prob.at(i);
  }
  return index.at(i);
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


arma::uvec sample_int(const unsigned & n,
                      const arma::uvec & index,
                      const arma::vec & prob) {
  arma::uvec ans(n);
  for (arma::uword k = 0; k < n; k++) {
    ans.at(k) = sample_int_1(index, prob);
  }
  return ans;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void merge_grids(arma::field<arma::vec> & g_loc, 
                 arma::field<arma::vec> & h_loc,
                 arma::vec & grid, 
                 arma::vec & h) {
  arma::uword n_grid = g_loc.n_rows;
  if (n_grid > 1) {
    grid.reset();
    h.reset();
    for (arma::uword k = 0; k < n_grid - 1; k++) {
      g_loc.at(k).resize(g_loc.at(k).n_rows - 1);
      grid = join_vert(grid, g_loc.at(k));
      h = join_vert(h, h_loc.at(k));
    }
    grid = join_vert(grid, g_loc.at(n_grid - 1));
    h = join_vert(h, h_loc.at(n_grid - 1));
  } else {
    grid = g_loc.at(0);
    h = h_loc.at(0);
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


arma::uvec n_non_zero_by_row(const arma::mat & map) {
  arma::uvec ans = arma::zeros<arma::uvec>(map.n_rows);
  arma::uvec indices;
  for (arma::uword k = 0; k < map.n_rows; k++){
    indices = find(map.row(k) > 0.);
    ans.at(k) = indices.n_rows;
  }
  return ans;
}

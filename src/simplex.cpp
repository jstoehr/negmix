// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

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
                   arma::vec & b) {
  // --- Dimension
  arma::uword k_p = w_p.n_rows;
  arma::uword k_n = w_n.n_rows;
  
  arma::uword i,j, n_pairs = a_star.n_rows;
  
  // --- Setting linear constraint
  arma::uword n_var = 2 * n_pairs;
  arma::uword n_eq = k_p + k_n + n_pairs;
  
  x = arma::zeros<arma::rowvec>(n_var);
  x.head(n_pairs).fill(1. - delta);
  x.tail(n_pairs).fill(-1.);
  
  // --- Bound
  b = join_vert(w_p, w_n, arma::zeros(n_pairs));
  
  // --- Constraint
  A = arma::zeros<arma::mat>(n_eq, n_var);
  
  for (arma::uword k = 0; k < n_pairs; k++) {
    i = pair.at(0, k);
    j = pair.at(1, k) + k_p;
    // --- Constraint on positive weights: the sum over the pair involving the
    // --- same positive component remain lower that the weight of the later
    A.at(i, k) = 1.;
    // --- Constraint on negative weights: the sum over the pair involving the
    // --- same positive component remain lower that the weight of the later
    A.at(j, n_pairs + k) = 1.;
    // --- Well defined pair
    A.at(k_p + k_n + k, k) = -1.;
    A.at(k_p + k_n + k, n_pairs + k) = a_star.at(k);
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void do_pivot(const arma::uword & i_pivot_row,
              const arma::uword & i_pivot_col,
              arma::rowvec & x,
              arma::mat & A,
              arma::vec & b,
              double & value) {
  
  // --- C++ version of the pivot function from boot package
  
  // --- Storing pivot quantities before update
  double pivot_value = A.at(i_pivot_row, i_pivot_col);
  double lbd_pivot = x.at(i_pivot_col);
  arma::vec pivot_col = A.col(i_pivot_col);
  double temp;
  
  // --- Update of the pivot row
  for (arma::uword j = 0; j < A.n_cols; j++) {
    if (A.at(i_pivot_row, j) != 0.) {
      A.at(i_pivot_row, j) /= (-1. * pivot_value);
    }
  }
  b.at(i_pivot_row) /= (-1. * pivot_value);
  
  A.at(i_pivot_row, i_pivot_col) = 1./pivot_value;
  
  // --- Update all rows except the pivot row 
  // --- Operation have been updated to handle
  // --- the update of A and b
  for (arma::uword i = 0; i < A.n_rows; i++) {
    if (i != i_pivot_row) {
      temp = A.at(i, i_pivot_col);
      if (temp != 0.) {
        for (arma::uword j = 0; j < A.n_cols; j++) {
          if (A.at(i_pivot_row, j) != 0.) {
            A.at(i,j) += temp * A.at(i_pivot_row, j);
          }
        }
        b.at(i) += temp * b.at(i_pivot_row);
        if (b.at(i) < 0.) {
          // --- Dealing with rounding issues
          b.at(i) = 0.;
        }
      }
    }
  }
  
  for (arma::uword j = 0; j < A.n_cols; j++) {
    if (A.at(i_pivot_row, j) != 0.) {
      x.at(j) += lbd_pivot * A.at(i_pivot_row, j);
    }
  }
  
  // --- New value of the objective function
  value += lbd_pivot * b.at(i_pivot_row);
  
  // ---
  for (arma::uword i = 0; i < A.n_rows; i++) {
    if (i != i_pivot_row) {
      A.at(i, i_pivot_col) = pivot_col.at(i)/pivot_value;
    }
  }
  x.at(i_pivot_col) = lbd_pivot/pivot_value;
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


void simplex(const double & tol,
             arma::rowvec & x,
             arma::mat & A, 
             arma::vec & b,
             arma::uvec & basic_var,
             double & value) {
  // --- Dimension
  arma::uword n_var = A.n_cols;
  arma::uword n_eq = A.n_rows;
  
  // --- For tracking the solution
  arma::uvec r_var = arma::linspace<arma::uvec>(0, n_var - 1, n_var);
  
  // --- Containers
  value = 0.;
  double min_of_r = x.at(0), min_of_c, count = 1.;
  bool stop, first_min;
  arma::uword j_temp = 0, i_temp, i_pivot_col, i_pivot_row;
  
  // --- First pivot col and checking optimality
  for (arma::uword j = 1; j < n_var; j++) {
    if (x.at(j) < min_of_r) {
      min_of_r = x.at(j);
      j_temp = j;
    }
  }
  
  if (min_of_r < -1. * tol) {
    stop = false;
  } else {
    stop = true;
  }
  
  while (!stop and count <= (n_var + n_eq)) {
    // --- Pivot column
    i_pivot_col = j_temp;
    
    first_min = true;
    for (arma::uword i = 0; i < n_eq; i++) {
      if (A.at(i, i_pivot_col) < -1. * tol) {
        double r = -b.at(i)/A.at(i, i_pivot_col);
        if (first_min) {
          first_min = false;
          min_of_c = r;
          i_pivot_row = i;
        } else if (r < min_of_c) {
          min_of_c = r;
          i_pivot_row = i;
        }
      }
    }
    if (first_min) {
      Rcpp::stop(" Error in cpp_simplex: no pivot row found!\n");
    }

    do_pivot(i_pivot_row, i_pivot_col, x, A, b, value);
    
    i_temp = basic_var.at(i_pivot_row);
    basic_var.at(i_pivot_row) = r_var.at(i_pivot_col);
    r_var.at(i_pivot_col) = i_temp;
    
    // --- First pivot col and checking optimality
    min_of_r = x.at(0);
    j_temp = 0;
    for (arma::uword j = 1; j < n_var; j++) {
      if (x.at(j) < min_of_r) {
        min_of_r = x.at(j);
        j_temp = j;
      }
    }
    // --- Checking if simplex is need
    if (min_of_r < -1. * tol) {
      stop = false;
    } else {
      stop = true;
    }
    ++count;
  }
  
  if (!stop) {
    Rcpp::stop(" Error in cpp_simplex: method did not converge!\n");
  }
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


Rcpp::List simplex_tab(const arma::vec & w_p,
                       const arma::vec & w_n,
                       const arma::umat & pair,
                       const arma::vec & a_star,
                       const double & delta) {
  arma::mat A;
  arma::vec b;
  arma::rowvec x;

  one_stage_tab(w_p, w_n, pair, a_star, delta,
                x, A, b);

  return Rcpp::List::create(
    Rcpp::Named("x", x),
    Rcpp::Named("A", A),
    Rcpp::Named("b", b)
  );
}


// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------


Rcpp::List rcpp_simplex(const arma::rowvec & x,
                        const arma::mat & A,
                        const arma::vec & b,
                        const arma::uvec & basic,
                        const double & tol) {
  arma::mat A_sol = -1. * A;
  arma::vec b_sol = b;
  arma::rowvec x_sol = x;
  arma::uvec basic_sol = basic - 1;

  double value;

  simplex(tol, x_sol, A_sol, b_sol, basic_sol, value);

  return Rcpp::List::create(
    Rcpp::Named("b_sol", b_sol),
    Rcpp::Named("A_sol", A_sol),
    Rcpp::Named("x_sol", x_sol),
    Rcpp::Named("basic_var", basic_sol),
    Rcpp::Named("value", value)
  );
}
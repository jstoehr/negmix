// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_zero.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

bool newton_raphson(f_zero & fun, 
                    double & x, 
                    double & f, 
                    double & df,
                    const int & maxit_0, 
                    const double & eps_0) {
  bool out = fun.grad(x, f, df);
  int count = 1;
  while (abs(f) > eps_0 & count < maxit_0) {
    x = x - f/df;
    out = fun.grad(x, f, df);
    if (out){
      break;
    }
    ++count;
  }
  if (count == maxit_0) {
    // printf("x %f f %f df %f\n", x, f, df);
    Rcpp::stop(" Newton Raphson method failed\n");
  }
  return out;
}

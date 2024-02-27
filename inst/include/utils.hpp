#ifndef INST_INCLUDE_UTILS_HPP_
#define INST_INCLUDE_UTILS_HPP_

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double inv_logit(const double & x, const double & inf, const double & sup,
                        const double & lambda) {
  return (sup - inf)/(1. + exp(-x/lambda)) + inf;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double logit(const double & x, const double & inf, const double & sup,
                    const double & lambda) {
  return -lambda * log((sup - inf)/(x - inf) - 1.);
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double grad_inv_logit(const double & x, const double & inf, const double & sup,
                             const double & lambda) {
  return (sup - inf)/(2. + exp(-x/lambda) + exp(x/lambda));
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double grad_transform_right(const double & x, const double & inf, const double & lambda) {
  return 1./lambda * exp(x/lambda);
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double log_sum_exp(const arma::vec & w, const arma::vec & x, arma::vec & w_exp) {
  double shift = x.max();
  arma::vec y = x - shift;
  w_exp = w % exp(y);
  return log(arma::sum(w_exp)) + shift;
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline void grid(const double & a,
                 const double & eps,
                 const double & b, 
                 arma::vec & g) {
  if (a + eps > b) {
    g = arma::linspace(a, b, 2);
  } else {
    g = arma::regspace(a, eps, b);
    g.at(g.n_rows - 1) = b;
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline void set_mono_comp(const arma::vec & base, 
                          const double & m, 
                          arma::ivec & mono) {
  mono.set_size(base.n_rows - 1);
  for (arma::uword k = 0; k < base.n_rows - 1; k++) {
    if (base.at(k + 1) <= m) {
      mono.at(k) = 1;
    } else {
      mono.at(k) = -1;
    }
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

arma::uword sample_int_1(const arma::uvec & index,
                         const arma::vec & prob);

arma::uvec sample_int(const unsigned & n,
                      const arma::uvec & index,
                      const arma::vec & prob);

void merge_grids(arma::field<arma::vec> & g_loc, 
                 arma::field<arma::vec> & h_loc,
                 arma::vec & grid, 
                 arma::vec & h);

arma::uvec n_non_zero_by_col(const arma::mat & map);

#endif /* INST_INCLUDE_UTILS_HPP_ */
// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include <f_model.hpp>
#include <simplex.hpp>


arma::mat f_model::map_pairs(const double & p_limit) {
  double a;
  double cond_a = 1./(1. - p_limit);
  arma::mat map = arma::zeros(k_p, k_n);
  
  for (arma::uword i = 0; i < k_p; i++) {
    for (arma::uword j = 0; j < k_n; j++) {
      if (this->is_valid_pair(i, j)) {
        a = this->a_star_pair(i, j);
        if(std::isinf(a) or (a > cond_a) or (a == 1.)) {
          map.at(i, j) = 0.;
        } else {
          map.at(i, j) = a;
        }
      }
    }
  }
  return map;
}


void f_model::build_pairs(const double & delta, const double & tol_simplex) {
  arma::mat map = this->map_pairs(delta);
  // --- Consider valid pairs only
  arma::uvec indices_star = find(map > 0.);
  arma::umat pair = arma::ind2sub(arma::size(map), indices_star);
  arma::vec a_star = map.elem(indices_star);
  
  if (a_star.n_rows == 0) {
    Rcpp::stop(" algorithm stopped: model does not contain possible pairs with mean acceptance probability lower than delta.\n");
    std::exit(0);
  }
  
  // --- Setting simplex table
  arma::mat A;
  arma::vec b;
  arma::rowvec x;
  
  one_stage_tab(w_p, w_n, pair, a_star, delta, x, A, b);
  
  // --- Simplex method
  arma::uword n_eq = A.n_rows;
  arma::uword n_var = A.n_cols;
  arma::uvec basic_var = arma::linspace<arma::uvec>(n_var, n_var + n_eq - 1, n_eq);;
  double value = 0.;
  A *= -1.;
  
  simplex(tol_simplex, x, A, b, basic_var, value);
  
  // --- Reading the result
  arma::vec sol = arma::zeros(n_var + n_eq);
  sol.elem(basic_var) = b;
  
  // --- Creating the pairs 
  arma::uword n_pairs = pair.n_cols;
  
  // --- Positive part of pairs
  w_pair_p = sol.subvec(0, n_pairs - 1);
  // --- Negative part of pairs
  w_pair_n = sol.subvec(n_pairs, n_var - 1);
  
  arma::uvec to_keep = find(w_pair_p > tol_simplex);
  w_pair_p = w_pair_p.elem(to_keep);
  w_pair_n = w_pair_n.elem(to_keep);
  list_pair = pair.cols(to_keep);
  
  // --- Setting the vector of proba to pick a pair
  p_accept_pair.set_size(to_keep.n_rows);
  // --- Initializing residuals weights
  r_p = w_p;
  r_n = w_n;
  
  arma::uword i,j;
  for (arma::uword k = 0; k < to_keep.n_rows; k++) {
    i = list_pair.at(0, k);
    j = list_pair.at(1, k);
    
    r_p.at(i) -= w_pair_p.at(k);
    if (r_p.at(i) < 0.) {
      // --- Numeric 0 issue
      r_p.at(i) = 0.;
    }
    
    r_n.at(j) -= w_pair_n.at(k);
    if (r_n.at(j) < 0.) {
      // --- Numeric 0 issue
      r_n.at(j) = 0.;
    }
    
    p_accept_pair.at(k) = 1. - w_pair_n.at(k)/w_pair_p.at(k);
  }
  // ---
  list_r_p = find(r_p != 0.);
  r_p = r_p.elem(list_r_p);
  
  list_r_n = find(r_n != 0.);
  r_n = r_n.elem(list_r_n);
  
  // --- Related outputs
  c_norm = w_pair_p - w_pair_n;
  this->set_pairs_param();
}


arma::mat f_model::augment_n_pair(const arma::uword & n_pair, const double & p_limit) {
  arma::uword k;
  arma::mat map = this->map_pairs(p_limit);
  arma::uvec count = n_non_zero_by_col(map);
  while (arma::sum(count) < n_pair) {
    // --- Replace a positive weight component by a positive weight component
    // --- that can balance more negative weight component
    k = count.index_min();
    this->set_rand_par_p(k);
    map = this->map_pairs(p_limit);
    count = n_non_zero_by_col(map);
  } 
  return map;
}


void f_model::ctrl_accept_rate(const double & a_ref, const double & p_max) {
  double coef = 1.;
  double temp = p_max * arma::sum(w_p) - 1.;
  if (temp < 0) {
    coef = arma::randu() * temp / (a_ref * (1. - p_max) - 1.);
  }
  double w_temp = coef / (a_ref - 1.);
  w_p.at(0) += w_temp * a_ref;
  w_n.at(0) += w_temp;
}


void f_model::create_benchmark_by_2(const arma::uword & k_init,
                                    const arma::uword & k_pair,
                                    const arma::uword & k_single,
                                    const arma::uword & n_pair,
                                    const double & p_min,
                                    const double & p_max,
                                    const double & p_star,
                                    const double & p_common,
                                    const double & p_rank,
                                    const double & p_limit,
                                    const int & maxit_0,
                                    const double & eps_0) {
  // printf(" Start rendering model by pair\n");
  // --- k_p positive weight components
  // --- k_n negative weight components
  // --- n_pair number of natural pairs
  // --- k_pair <= k_p number of positive component involved in natural pairing
  // --- k_init <= k_pair and k_n positive weight components used to init negative parameters
  // --- k_single <= k_p - k_pair positive component that cannot balance any negative weight components
  
  arma::uword k_temp = k_n + k_pair - k_init, i, j;
  arma::umat l_init(2, k_temp);
  arma::uvec init_par(k_temp);
  arma::vec a_pair(n_pair);
  
  for (arma::uword k = 0; k < k_n; k++) {
    // --- Negative weight component based on available positive weight component
    if (k < k_init) {
      l_init.at(0, k) = k;
      l_init.at(1, k) = k;
    } else {
      l_init.at(0, k) = arma::randi<arma::uword>(arma::distr_param(0, k_init - 1));
      l_init.at(1, k) = k;
    }
    init_par.at(k) = 1;
  }
  
  for (arma::uword k = k_init; k < k_pair; k++) {
    // --- Positive weight component based on available negative weight component
    l_init.at(0, k + k_n - k_init) = k;
    l_init.at(1, k + k_n - k_init) = arma::randi<arma::uword>(arma::distr_param(0, k_n - 1));
    init_par.at(k + k_n - k_init) = 2;
  }
  
  // --- First pair is set to control the overall acceptance rate
  // printf(" Setting the reference pair...");
  double p = arma::randu(arma::distr_param(p_min, p_max));
  bool is_star = (arma::randu() < p_star);
  bool has_common = (arma::randu() < p_common);
  this->render_pair_from_pc(0, 0, p, is_star, has_common, maxit_0, eps_0);
  double a_ref = 1./(1. - p);
  // printf(" DONE!\n");
  
  // --- Initialization of negative and pair parameters
  // printf(" Initialization components paramater...");
  for (arma::uword k = 1; k < k_temp; k++) {
    i = l_init.at(0, k);
    j = l_init.at(1, k);
    is_star = (arma::randu() < p_star);
    has_common = (arma::randu() < p_common);
    if (arma::randu() < p_rank) {
      p = arma::randu(arma::distr_param(p_min, p_max));
    } else {
      p = arma::randu();
    }
    if (init_par.at(k) == 1) {
      this->render_pair_from_pc(i, j, p, is_star, has_common, maxit_0, eps_0);
    } else if (init_par.at(k) == 2) {
      this->render_pair_from_nc(i, j, p, is_star, has_common, maxit_0, eps_0);
    }
    a_pair.at(k) = 1./(1. - p);
  }
  // printf(" DONE!\n");
  
  // --- Correcting parameters to get exactly n_pairs
  // printf(" Correcting number of pairs...");
  if (k_temp < n_pair) {
    arma::mat map = this->augment_n_pair(n_pair, p_limit);
    for (arma::uword k = 1; k < k_temp; k++) {
      i = l_init.at(0, k);
      j = l_init.at(1, k);
      map.at(i, j) = 0.;
    }
    arma::uvec indices = find(map > 0.);
    arma::uvec new_index = sample_int(n_pair - k_temp, 
                                      indices, arma::ones(indices.n_rows)/indices.n_rows);
    // --- Update list of pairs
    l_init = arma::join_horiz(l_init, arma::ind2sub(arma::size(map), new_index));
    // --- Compute a for new pairs
    arma::vec a_star = map.elem(new_index);
    for (arma::uword k = k_temp; k < n_pair; k++) {
      if(arma::randu() < p_star) {
        a_pair.at(k) = a_star.at(k - k_temp);
      } else {
        double a_temp;
        if (a_star.at(k - k_temp) < 1./(1. - p_max) and arma::randu() < p_rank) {
          a_temp = 1./(1. - p_max);
        } else {
          a_temp = 1./(1. - p_limit);
        }
        a_pair.at(k) = arma::randu(arma::distr_param(a_star.at(k - k_temp), a_temp));
      }
    }
  }
  // printf(" DONE!\n");
  
  // --- Computing weight for a mixture (minus the first pair 
  // --- used to control the acceptance rate)
  arma::vec w_mix = arma::randu(n_pair);
  w_mix.at(0) = 0.;
  w_mix /= arma::sum(w_mix);
  double w_temp;
  
  for (arma::uword k = 1; k < n_pair; k++) {
    i = l_init.at(0, k);
    j = l_init.at(1, k);
    w_temp = w_mix.at(k)/(a_pair.at(k) - 1.);
    w_p.at(i) += w_temp * a_pair.at(k);
    w_n.at(j) += w_temp;
  }
  
  this->ctrl_accept_rate(a_ref, p_max);
  this->normalize_weight();
  
  // --- Dealing with single positive components
  if (k_p - k_pair > 0) {
    // printf(" Adding single positive weight components...");
    arma::uword k_p_1 = k_p - k_pair - k_single;
    this->add_single_comp(k_p_1, k_single);
    double coef = arma::randu() * (p_max * sum(w_p) - 1.) / (1. - p_max);
    arma::vec prop = arma::randu(k_p - k_pair);
    prop /= arma::sum(prop);
    w_p.tail(k_p - k_pair) = coef * prop;
    this->normalize_weight();
    // printf(" DONE!\n");
  }
  
  this->set_alt_param();
}


void f_model::create_benchmark(const arma::uword & k_init,
                               const double & p_min,
                               const double & p_max,
                               const double & p_star,
                               const double & p_common,
                               const double & p_rank,
                               const int & maxit_0, 
                               const double & eps_0, 
                               const int & maxit, 
                               const double & eps_f, 
                               const double & eps_g) {
  // printf(" Start rendering model\n");
  // --- Containers
  arma::uvec l_pc(k_n);
  arma::uword i;
  arma::vec a_s(k_n), a(k_n);
  
  for (arma::uword k = 0; k < k_n; k++) {
    // --- Negative weight component based on available positive weight component
    if (k < k_init) {
      l_pc.at(k) = k;
    } else {
      l_pc.at(k) = arma::randi<arma::uword>(arma::distr_param(0, k_init - 1));
    }
  }
  
  // --- First pair is set to control the overall acceptance rate
  // printf(" Setting the reference pair...");
  double p = arma::randu(arma::distr_param(p_min, p_max));
  bool is_star = (arma::randu() < p_star);
  bool has_common = (arma::randu() < p_common);
  this->render_pair_from_pc(0, 0, p, is_star, has_common, maxit_0, eps_0);
  double a_ref = 1./(1. - p);
  a_s.at(0) = this->a_star_pair(0, 0);
  // printf(" DONE!\n");
  
  // --- Initializing negative weight components parameters
  // printf(" Initialization negative weight components paramater...");
  for (arma::uword k = 1; k < k_n; k++) {
    i = l_pc.at(k);
    has_common = (arma::randu() < p_common);
    if (arma::randu() < p_rank) {
      p = arma::randu(arma::distr_param(p_min, p_max));
    } else {
      p = arma::randu();
    }
    this->render_pair_from_pc(i, k, p, true, has_common, maxit_0, eps_0);
    a_s.at(k) = this->a_star_pair(i, k);
  }
  // printf(" DONE!\n");
  
  // --- Compute weight to balance negetive part
  // printf(" Computing weights...");
  double w_n_ref = median(1. / (a_s - 1.));
  arma::vec coef = arma::randu(k_n);
  coef /= arma::sum(coef);
  arma::mat w_p_temp(k_p, k_n);
  for (arma::uword k = 1; k < k_n; k++) {
    i = l_pc.at(k);
    w_n.at(k) = coef.at(k) * w_n_ref;
    w_p_temp.at(i, k) = w_n.at(k) * arma::randu() * a_s.at(k);
    w_p_temp.col(k) = this->w_valid_mixt(i, k, w_p_temp.at(i, k), w_n.at(k), 
                 maxit_0, eps_0, maxit, eps_f, eps_g);
  }
  // printf(" DONE!\n");
  
  w_p = arma::sum(w_p_temp, 1);
  this->normalize_weight();
  
  // --- Add first component with proper coefficient to control acceptance rate
  this->ctrl_accept_rate(a_ref, p_max);
  this->normalize_weight();
  this->set_alt_param();
}

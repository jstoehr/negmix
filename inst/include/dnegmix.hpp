#ifndef INST_INCLUDE_DNEGMIX_HPP_
#define INST_INCLUDE_DNEGMIX_HPP_

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double log_gampdf_1(const double & x,
                           const double & alpha,
                           const double & beta) {
  if (alpha == 1.) {
    return log(beta) - beta * x;
  } else {
    return (alpha - 1.) * log(x) + alpha * log(beta) - beta * x - std::lgamma(alpha);
  }
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

arma::vec dnormix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & mean, 
                  const arma::vec & sd);

double dgammix_1(const double & x,
                 const arma::vec & w,
                 const arma::vec & alpha,
                 const arma::vec & beta);

arma::vec dgammix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & alpha, 
                  const arma::vec & beta); 

#endif /* INST_INCLUDE_DNEGMIX_HPP_ */
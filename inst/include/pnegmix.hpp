#ifndef INST_INCLUDE_PNEGMIX_HPP_
#define INST_INCLUDE_PNEGMIX_HPP_

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

arma::vec pnormix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & mean, 
                  const arma::vec & sd);
  
arma::vec pgammix(const arma::vec & x, 
                  const arma::vec & w, 
                  const arma::vec & alpha, 
                  const arma::vec & inv_beta);

#endif /* INST_INCLUDE_PNEGMIX_HPP_ */
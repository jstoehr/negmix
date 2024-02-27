#ifndef INST_INCLUDE_QNEGMIX_HPP_
#define INST_INCLUDE_QNEGMIX_HPP_

#include <f_model.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

void inv_cdf_method(f_model & f,
                    arma::vec & ans,
                    double & cpu_time,
                    const double & precision = 1e-10,
                    const unsigned & n_bins = 500);

#endif /* INST_INCLUDE_QNEGMIX_HPP_ */
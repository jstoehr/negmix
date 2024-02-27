#ifndef INST_INCLUDE_A_STAR_GAM_HPP_
#define INST_INCLUDE_A_STAR_GAM_HPP_

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double x_star_norm(const double & mean_p,
                          const double & sd_p,
                          const double & mean_n,
                          const double & sd_n) {
  return (sd_p * sd_p * mean_n - sd_n * sd_n * mean_p) / (sd_p * sd_p - sd_n * sd_n);
}

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

inline double x_star_gam(const double & alpha_p,
                         const double & beta_p, 
                         const double & alpha_n, 
                         const double & beta_n) {
  return (alpha_p - alpha_n)/(beta_p - beta_n);
};

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

double a_star_norm(const double & mean_p, 
                   const double & sd_p, 
                   const double & mean_n, 
                   const double & sd_n,
                   const bool & do_log);

double a_star_gam(const double & alpha_p, 
                  const double & beta_p, 
                  const double & alpha_n, 
                  const double & beta_n,
                  const bool & do_log);

#endif /* INST_INCLUDE_A_STAR_GAM_HPP_ */
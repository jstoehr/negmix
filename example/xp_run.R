rm(list = ls())
library(negmix)

case <- "norm"
if (case == "norm") {
  family <- "normal"
} else if (case == "gam") {
  family <- "gamma"
}
dir <- paste0("negmix_", case, "/")
dir_res <- paste0(dir, "result/")
if (!dir.exists(dir_res)) {
  dir.create(dir_res, recursive = TRUE)
}

# --- Model setting
set_config <- 1:4
set_param <- 1:2
xp <- 1:50
set_p <- c(0., 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 3e-1)

# --- Experiment setting
n_sample <- c(1e1, 1e2, 1e3, 1e4)
tol <- 1e-10
delta <- c(.4, .6, .8)
set_eps <- c(.1, .2, .5, 1.)
n_eps <- 0
for (d in delta){
  n_eps <- n_eps + sum(set_eps < (1 - d)/d)
}

# --- Running experiment
n_exp <- length(set_config) * length(set_param) * length(xp) * length(set_p[-1])
(n_exp <- n_exp * length(n_sample) * n_eps)
pb <- txtProgressBar(min = 0, max = n_exp, 
                     style = 3, width = 50, char = "=")

i <- 0
overwrite_ar <- FALSE
overwrite_inv <- FALSE
overwrite_xp <- FALSE
for (n in n_sample) {
  for (param in set_param) {
    for (setting in set_config) {
      for (j in 2:length(set_p)) {
        p_min <- set_p[j - 1]
        p_max <- set_p[j]
        for (seed in xp) {
          name_config <- paste0(
            "rand_neg_", case, "_set", setting, "_param", param, 
            "_seed", seed, "_pmin", p_min, "_pmax", p_max)
          
          filename <- paste0(dir, name_config, ".RData")
          filename_ar <- paste0(
            dir_res, name_config, "_sample_n", n,
            "_ar_vanilla.RData"
          )
          filename_inv <- paste0(
            dir_res, name_config, "_sample_n", n,
            "_cdf_inv.RData"
          )
          
          # --- Loading model
          if (file.exists(filename)) {
            load(filename)
            if (!file.exists(filename_ar) || overwrite_ar) {
              out_vanilla <- rnegmix(
                n, par, family = family, delta = NULL,
                method = "vanilla"
              )
              save(out_vanilla, file = filename_ar)
              overwrite_result <- TRUE
            }
            
            if (!file.exists(filename_inv) || overwrite_inv) {
              out_inv <- rnegmix(
                n, par, family = family, delta = NULL,
                method = "inv_cdf", control = list(breaks = 500)
              )
              save(out_inv, file = filename_inv)
              overwrite_result <- TRUE
            }
            
            for (d in delta) {
              eps <- set_eps[set_eps < (1 - d)/d]
              for (e in eps) {
                # --- File where to save outputs
                filename_xp <- paste0(
                  dir_res, name_config, "_sample_n", n,
                  "_d", d, "_eps", e, ".RData"
                )
                
                if (!file.exists(filename_xp) || overwrite_xp) {
                  rnegmix(
                    n, par, family = family, delta = d,
                    method = "stratified", 
                    control = list(eps_d = e, optim = FALSE, use_mono = TRUE, n_points = 100)
                  )
                  save(out_strat, file = filename_xp)
                  overwrite_result <- TRUE
                }
                # Sets the progress bar to the current state
                i <- i + 1
                setTxtProgressBar(pb, i)
              }
              # --- End of experiment on epsilon
            }
            # --- End of experiment on delta
            # --- End of all experiments for that model
          } else {
            print("No file: ", name_config, "\n")
          }
        }
      }
    }
  }
}
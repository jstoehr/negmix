library(negmix)

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


overwrite_result <- FALSE
pb <- txtProgressBar(min = 0, max = n_exp,
                     style = 3, width = 50, char = "=")
i <- 0
output <- c()
r_n <- c()
for (n in n_sample) {
  # --- Merging the result
  result_name <- paste0(dir_res, "rand_neg_", case, "_sample_n", n, "_result.RData")
  if (overwrite_result || !file.exists(result_name)) {
    result <- c()
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
            load(filename)
            load(filename_ar)
            load(filename_inv)
            
            for (d in delta) {
              eps <- set_eps[set_eps < (1 - d)/d]
              map <- cpp_map_pairs(d, par, family = family)
              for (e in eps) {
                filename_xp <- paste0(
                  dir_res, name_config, "_sample_n", n,
                  "_d", d, "_eps", e, ".RData"
                )
                load(filename_xp)
                
                res_config <- data.frame(
                  n = n,
                  set = setting,
                  param = param,
                  xp = seed,
                  p_min = p_min,
                  p_max = p_max,
                  cat = max(which(set_p < 1/sum(par$w_p))),
                  k_p = length(par$w_p),
                  k_n = length(par$w_n),
                  prob_v = out_vanilla$prob_accept,
                  freq_v = out_vanilla$freq_accept,
                  time_ar = out_vanilla$time,
                  time_inv = out_inv$time,
                  delta = d,
                  eps = e,
                  n_pairs = ifelse(sum(map) > 0,
                                   ncol(out_strat$pairs), NA),
                  n_max = sum(map > 0),
                  prob_s = ifelse(sum(map) > 0,
                                  out_strat$prob_accept, NA),
                  freq_s = ifelse(sum(map) > 0,
                                  out_strat$freq_accept, NA),
                  sum_r_p = ifelse(length(out_strat$r_p) > 0,
                                   sum(out_strat$r_p), 0),
                  sum_r_n = ifelse(length(out_strat$r_n) > 0,
                                   sum(out_strat$r_n), 0),
                  time = ifelse(sum(map) > 0,
                                out_strat$time, NA),
                  precision = 1e-10
                )
                
                result <- rbind(result, res_config)
                # Sets the progress bar to the current state
                i <- i + 1
                setTxtProgressBar(pb, i)
              }
              # --- End of experiment on epsilon
            }
            # --- End of experiment on delta
            # --- End of all experiments for that model
          }
        }
      }
    }
    save(result, file = result_name)
    output <- rbind(output, result)
  } else {
    load(result_name)
    output <- rbind(output, result)
  }
}


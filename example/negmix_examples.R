rm(list = ls())
library(negmix)

case <- "gam"
dir <- paste0("example/negmix_", case, "/")
if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

# ------------- EXAMPLES WITH NATURAL PAIRING ------------

for (ex in 1:2) {
  for (inflated in c(TRUE, FALSE)) {
    filename <- paste0(dir, "negmix_", case, "_ex", ex, "_inf", inflated * 1., ".RData")
    if (ex == 1 & case == "gam") {
      # --- Exemple n°1
      shape_p <- seq(1, 5, .1)
      n_p <- length(shape_p)
      rate_p <- seq(0.25, 2, length.out = n_p)
      shape_n <- shape_p + .01
      rate_n <- rate_p + .01
    } else if (ex == 2 & case == "gam") {
      # --- Exemple n°2
      shape_p <- rep(1, 10)
      n_p <- length(shape_p)
      rate_p <- seq(0.5, 4, length.out = n_p)
      shape_n <- shape_p
      rate_n <- rate_p + runif(n_p)
    } else if (ex == 1 & case == "norm") {
      # --- Exemple n°1
      mean_p <- seq(0, 10, .2)
      n_p <- length(mean_p)
      sd_p <- seq(.25, 1, length.out = n_p)
      mean_n <- mean_p + .01
      sd_n <- sd_p - .01
      
    } else if (ex == 2 & case == "norm") {
      # --- Exemple n°2
      mean_p <- c(-5, -4.8, -4.6, 0, 4.8, 5)
      n_p <- length(mean_p)
      sd_p <- rep(1, n_p)
      mean_n <- c(-4.9, -4.9, -4.7, 0, 4.9, 4.9)
      sd_n <- rep(.9, n_p)
    }
    
    map_pairs_m <- diag(1, n_p, n_p)
    a_star <- rep(0, n_p)
    if (case == "gam") {
      par <- list(w_p = NULL, shape_p = shape_p, rate_p = rate_p,
                  w_n = NULL, shape_n = shape_n, rate_n = rate_n)
      for (i in seq_len(n_p)) {
        a_star[i] <- a_star_gam(shape_p[i], rate_p[i], shape_n[i], rate_n[i], FALSE)
      }
    } else if (case == "norm") {
      par <- list(w_p = NULL, mean_p = mean_p, sd_p = sd_p,
                  w_n = NULL, mean_n = mean_n, sd_n = sd_n)
      for (i in seq_len(n_p)) {
        a_star[i] <- a_star_norm(mean_p[i], sd_p[i], mean_n[i], sd_n[i], FALSE)
      }
    }
    # --- Computing weights
    par$w_p <- a_star/(a_star - 1)
    par$w_n <- 1/(a_star - 1)
    if (inflated) {
      par$w_n <- par$w_n * par$w_p
      par$w_p <- par$w_p * par$w_p
    }
    c_norm <- sum(par$w_p) - sum(par$w_n)
    par$w_p <- par$w_p/c_norm
    par$w_n <- par$w_n/c_norm
    save(par, map_pairs_m, file = filename)
  }
}

# ------------- RANDOM EXAMPLES ------------

overwrite <- FALSE
set_config <- 1:4
set_param <- 1:2
xp <- 50
set_p <- c(0., 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 3e-1)

n_exp <- length(set_config) * length(set_param) * xp * length(set_p[-1])
pb <- txtProgressBar(min = 0, max = n_exp, 
                     style = 3, width = 50, char = "=")
i <- 0
for (setting in set_config) {
  for (param in set_param) {
    
    file_set <- paste0(dir, "rand_neg_", case, "_set", setting, "_param", param)
    out_set <- paste0(file_set, ".RData")
    
    if (setting == 1) {
      k_p_min <- 5
      k_p_max <- 10
      k_n_max <- 2
    } else if (setting == 2) {
      k_p_min <- 10
      k_p_max <- 30
      k_n_max <- 2
    } else if (setting == 3) {
      k_p_min <- 30
      k_p_max <- 50
      k_n_max <- 2
    } else if (setting == 4) {
      k_p_min <- 50
      k_p_max <- 100
      k_n_max <- 2
    } 
    
    for (j in 2:length(set_p)) {
      p_min <- set_p[j - 1]
      p_max <- set_p[j]
      for (seed in seq_len(xp)) {
        filename <- paste0(
          file_set, "_seed", seed, 
          "_pmin", p_min, "_pmax", p_max, ".RData"
        )
        if (!file.exists(filename) || overwrite) {
          set.seed(seed)
          k_p <- sample(k_p_min:k_p_max, 1)
          k_n <- sum(sample(1:k_n_max, k_p, TRUE))
          k_p_init <- ifelse(param == 1, sample(1:min(k_n, k_p), 1), k_p)
          if (param == 1) {
            k_p_pair <- ifelse(k_p_init < k_p, sample(k_p_init:k_p, 1), k_p)
            n_pair <- k_n + k_p_pair - k_p_init
            n_single <- sample(0:(k_p - k_p_pair), 1)
            control <- list(
              p_min = p_min,
              p_max = p_max,
              p_star = 1.,
              k_p_init = k_p_init,
              k_p_pair = k_p_pair,
              n_pair = n_pair,
              n_single = n_single
            )
          } else {
            control <- list(
              p_min = p_min,
              p_max = p_max,
              p_star = 1.,
              k_p_init = k_p_init
            )
          }
          
          par <- render_negmix(
            k_p, k_n,
            family = case,
            by_pair = (param == 1),
            control = control
          )
          save(par, file = filename)
          # Sets the progress bar to the current state
          i <- i + 1
          setTxtProgressBar(pb, i)
        }
      }
    }
  }
}

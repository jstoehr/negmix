render_negmix <- function(k_p, k_n, family = c("normal", "gamma"), par = list(),
                          by_pair = FALSE, control = list()) {
  # --- Family check
  family <- match.arg(family)
  # --- Par check
  if (length(par) > 0) {
    names_par <- names(par)
    if (family == "normal") {
      names_model <- c("mean_p", "sd_p")
    } else if (family == "gamma") {
      names_model <- c("shape_p", "rate_p")
    }
    if (length(unknown <- names_par[!names_par %in% names_model])) {
      warning("unknown names in par: ", paste(unknown, collapse = ", "))
    }
    if (length(missing <- names_model[!names_model %in% names_par])) {
      stop("missing argument with no default value in par: ", paste(missing, collapse = ", "))
    }
    par <- par[names_model]
    if(!all(lengths(par) == length(par[[1]])) & min(lengths(par)) > 0) {
      stop("Expecting arguments of the same length in par")
    } else if(length(par[[1]]) > k_p) {
      stop("Expecting arguments of length lower than k_p [=", k_p, "]")
    } else if (is.null(control$k_init)) {
      control$k_init <- length(par[[1]])
    }
  }
  # --- Control check
  ctrl <- list(
    k_init = k_n,
    k_pair = k_n,
    k_single = 0,
    n_pair = k_n,
    p_min = 0.,
    p_max = 1.,
    p_star = 0.8,
    p_common = 0.2,
    p_rank = 1.,
    p_limit = .95,
    lambda = 1.,
    maxit_0 = 1e6,
    eps_0 = 1e-10,
    maxit = 300,
    eps_f = 1e-5,
    eps_g = 1e-6
  )
  names_ctrl <- names(ctrl)
  names_control <- names(control)
  ctrl[names_control] <- control
  if (length(unknown <- names_control[!names_control %in% names_ctrl])) {
    warning("unknown names in control: ", paste(unknown, collapse = ", "))
  }
  
  if (is.null(control$k_init)) {
    ctrl$k_init <- ifelse(by_pair, sample(1:k_n, 1), sample(1:k_p, 1))
  }
  if (is.null(control$k_pair)) {
    ctrl$k_pair <- ifelse(ctrl$k_init < k_p, sample(ctrl$k_init:k_p, 1), k_p)
  }
  if (is.null(control$n_pair)) {
    ctrl$n_pair <- k_n + ctrl$k_pair - ctrl$k_init
  } else {
    ctrl$n_pair - k_n > ctrl$k_pair - ctrl$k_init
  }
  
  
  if (by_pair) {
    if (ctrl$k_init > k_n) {
      warning("length of initial parameters [=", ctrl$k_init, "] is larger than k_n [=", k_n, "]. Reset to k_n")
      ctrl$k_init <- k_n
    }
    if (ctrl$k_pair < ctrl$k_init || ctrl$k_pair > k_p) {
      stop("invalid value for control$k_pair [=", ctrl$k_pair, "]. 
           Expecting value between control$k_init [=", ctrl$k_init,"] and k_p [=", k_p, "]")
    }
    if (ctrl$k_pair + ctrl$k_single > k_p) {
      stop("invalid value for control$k_single [=", ctrl$k_single, "]. 
           Expecting value between 0 and k_p - control$k_pair [=", k_p - ctrl$k_pair, "]")
    }
    if (ctrl$n_pair < k_n + ctrl$k_pair - ctrl$k_init) {
      stop("invalid value for control$n_pair [=", ctrl$n_pair, "]. 
           Cannot be lower than k_n + control$k_pair - control$k_init  [=", k_n + ctrl$k_pair - ctrl$k_init, "]")
    } else if (ctrl$n_pair > k_p * k_n) {
      stop("invalid value for control$n_pair. Cannot exceed k_p * k_n [=", k_p * k_n, "]")
    }
  } else if (length(ignored <- names_control[names_control %in% c("k_pair", "k_single", "n_pair")])) {
    warning("argument ignored in control when by_pair = FALSE: ", paste(ignored, collapse = ", "))
  }
  
  # --- par check
  k <- ifelse(by_pair, ctrl$k_init, k_p)
  if (family == "normal") {
    param <- list(mean_p = runif(k, 0., 20.), sd_p = rgamma(k, 3., 2.5))
  } else if (family == "gamma") {
    param <- list(shape_p = rgamma(k, 4., 0.5), rate_p = rgamma(k, 2., 0.7))
    if (runif(1L) < .5) {
      # --- 50% to have a model with either only shape_p < 1 or shape_p >= 1
      if (sum(param$shape_p >= 1.) == k_p) {
        a <- sample(param$shape_p, 1L)
        param$shape_p <- param$shape_p / (a + .1)
      } else if (sum(param$shape_p < 1.) > 0) {
        param$shape_p <- param$shape_p + 1. - min(param$shape_p) + runif(1L)
        param$shape_p <- sort(param$shape_p)
        mode <- seq(0.5, runif(1L, 3., 7.), length.out = k_p)
        param$rate_p <- (param$shape_p - 1.)/mode
      }
    }
  }
  
  if (length(par) == 0) {
    par <- param
  } else {
    n_ <- length(par[[1]])
    if ((by_pair & n_ < ctrl$k_init) || (!by_pair & n_ < k_p)) {
      cat("je suis ici\n")
      for (i in 1:length(param)) {
        param[[i]][1:n_] <- par[[i]]
      }
      par <- param
    } else if (by_pair & n_ > ctrl$k_init) {
      warning("length of initial parameters [=", ctrl$k_init, "] is larger than control$k_init [=", ctrl$k_init, "]. 
              Keeping only the first control$k_init elements of each vector in argument par")
      for (i in 1:n_) {
        par[[i]] <- par[[i]][1:ctrl$k_init]
      }
    }
  }
  
  return(lapply(render_model(k_p, k_n, family, par, by_pair, ctrl), as.numeric))
}

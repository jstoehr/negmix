render_negmix <- function(k_p, k_n, family = c("normal", "gamma"),
                          by_pair = FALSE, control = list()) {
  # --- Family check
  family <- match.arg(family)
  # --- Control check
  ctrl <- list(
    k_p_init = ifelse(by_pair, sample(1:min(k_n, k_p), 1), sample(1:k_p, 1)),
    k_p_pair = k_p,
    n_single = 0,
    n_pair = k_n,
    p_min = 0.,
    p_max = 0.95,
    p_star = 0.8,
    p_common = 0.2,
    p_rank = 1.,
    p_limit = .99,
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
  
  if (ctrl$k_p_init > k_p) {
    stop("invalid value for control$k_p_init [=", ctrl$k_p_init, "]. 
           Expecting value lower than k_p [=", k_p,"]")
  }
  if (by_pair) {
    if (is.null(control$k_p_pair)) {
      ctrl$k_p_pair <- ifelse(ctrl$k_p_init < k_p, sample(ctrl$k_p_init:k_p, 1), k_p)
    } else if (is.null(control$k_p_init)) {
      ctrl$k_p_init <- sample(1:ctrl$k_p_pair, 1)
    }
    if (is.null(control$n_pair)) {
      k_ <- ctrl$k_p_pair - ctrl$k_p_init + k_n
      ctrl$n_pair <- ifelse(
        k_ == k_n * (ctrl$k_p_pair - ctrl$k_p_init + 1), 
        k_, 
        sample(k_:(k_n * (ctrl$k_p_pair - ctrl$k_p_init + 1)), 1)
      )
    } else {
      if (is.null(control$k_p_pair) & is.null(control$k_p_init)) {
        if (ctrl$n_pair > k_n) {
          k_ <- ceiling(ctrl$n_pair/k_n - 1)
          diff <- ifelse(
            k_ >= k_p - 1, 
            k_p - 1, 
            sample(k_:(k_p - 1), 1)
          )
          ctrl$k_p_init <- ifelse(
            k_p - diff == 1, 
            1, 
            sample(1:(k_p - diff), 1)
          )
          ctrl$k_p_pair <- ctrl$k_p_init + diff
        }
      } else if (is.null(control$k_p_pair)) {
        k_ <- ceiling(ctrl$n_pair/k_n - 1 + ctrl$k_p_init)
        k_max <- ctrl$n_pair + ctrl$k_init - k_n
        ctrl$k_p_pair <- ifelse(
          k_ == k_max,
          k_max,
          sample(k_:k_max, 1)
        )
      } else if (is.null(control$k_p_init)) {
        k_ <- floor(ctrl$k_p_pair - ctrl$n_pair/k_n + 1)
        k_min <- k_n + ctrl$k_p_pair - ctrl$n_pair
        ctrl$k_p_init <- ifelse(
          k_ == k_min,
          k_min,
          sample(k_min:k_, 1)
        )
      }
    }
    if (is.null(control$n_single)) {
      ctrl$n_single <- sample(0:(k_p - ctrl$k_p_pair), 1)
    }
    if (ctrl$k_p_pair < ctrl$k_p_init || ctrl$k_p_pair > k_p) {
      stop("invalid value for control$k_p_pair [=", ctrl$k_p_pair, "]. 
           Expecting value between control$k_p_init [=", ctrl$k_p_init,"] and k_p [=", k_p, "]")
    }
    if (ctrl$k_p_pair + ctrl$n_single > k_p) {
      stop("invalid value for control$n_single [=", ctrl$n_single, "]. 
           Expecting value between 0 and k_p - control$k_p_pair [=", k_p - ctrl$k_p_pair, "]")
    }
    if (ctrl$n_pair < k_n + ctrl$k_p_pair - ctrl$k_p_init) {
      stop("invalid value for control$n_pair [=", ctrl$n_pair, "]. 
           Cannot include less than k_n + control$k_p_pair - control$k_p_init  [=", 
           k_n + ctrl$k_p_pair - ctrl$k_p_init, 
           "] two-component signed mixtures")
    } else if (ctrl$n_pair > k_n * (ctrl$k_p_pair - ctrl$k_p_init + 1)) {
      warning("cannot ensure to include more than 
              k_n * (control$k_p_pair - control$k_p_init + 1) [=", 
              k_n * (ctrl$k_p_pair - ctrl$k_p_init + 1), "] two-component signed mixtures")
    }
  } else if (length(ignored <- names_control[names_control %in% c("k_p_pair", "n_single", "n_pair")])) {
    warning("argument ignored in control when by_pair = FALSE: ", paste(ignored, collapse = ", "))
  }
  if (ctrl$p_max == 1.) {
    stop("invalid value for control$p_max [=", ctrl$p_max, "]. Expecting a value strictly lower than 1")
  }
  
  # --- 
  k <- ifelse(by_pair, ctrl$k_p_init, k_p)
  if (family == "normal") {
    param <- list(mean_p = runif(k, 0., 20.), sd_p = rgamma(k, 3., 2.5))
  } else if (family == "gamma") {
    param <- list(shape_p = rgamma(k, 4., 0.5), rate_p = rgamma(k, 2., 0.7))
    if (runif(1L) < .5) {
      # --- 50% to have a model with either only shape_p < 1 or shape_p >= 1
      if (sum(param$shape_p >= 1.) == k) {
        a <- sample(param$shape_p, 1L)
        param$shape_p <- param$shape_p / (a + .1)
      } else if (sum(param$shape_p < 1.) > 0) {
        param$shape_p <- param$shape_p + 1. - min(param$shape_p) + runif(1L)
        param$shape_p <- sort(param$shape_p)
        mode <- seq(0.5, runif(1L, 3., 7.), length.out = k)
        param$rate_p <- (param$shape_p - 1.)/mode
      }
    }
  }
  return(lapply(render_model(k_p, k_n, family, param, by_pair, ctrl), as.numeric))
}

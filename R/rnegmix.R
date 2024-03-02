rnegmix <- function(n, par, family = c("normal", "gamma"), delta = NULL,
                    method = c("stratified", "vanilla", "inv_cdf"), control = list()) {
  # --- Family check
  family <- match.arg(family)
  # --- Par check
  names_par <- names(par)
  if (family == "normal") {
    names_p <- c("w_p", "mean_p", "sd_p")
    names_n <- c("w_n", "mean_n", "sd_n")
  } else if (family == "gamma") {
    names_p <- c("w_p", "shape_p", "rate_p")
    names_n <- c("w_n", "shape_n", "rate_n")
  }
  names_model <- c(names_p, names_n)
  if (length(unknown <- names_par[!names_par %in% names_model])) {
    warning("unknown names in par: ", paste(unknown, collapse = ", "))
  }
  if (length(missing <- names_model[!names_model %in% names_par])) {
    stop("missing argument with no default value in par: ", paste(missing, collapse = ", "))
  }
  if (sum(par$w_n < 0) > 0) {
    stop("expecting positive values for par$w_n")
  } 
  if (!all(lengths(par[names_par %in% names_p]) == length(par$w_p))) {
    stop("expecting positive weight components parameters of same size in par")
  }
  if (!all(lengths(par[names_par %in% names_n]) == length(par$w_n))) {
    stop("expecting negative weight components parameters of same size in par")
  }
  # --- Method check
  method <- match.arg(method)
  if(is.null(delta) & method == "stratified") {
    stop("delta is missing with no default for the stratified method")
  } else if (!is.null(delta)) {
    if (method != "stratified") {
      warning("delta ignored for vanilla and inv_cdf methods")
    } else if (delta < 1/sum(par$w_p)) {
      warning("vanilla method performs better. Increase delta value")
    }
  }
  # --- Control check
  ctrl <- list(delta = delta,
               eps_d = 0.9 * (1. - delta)/delta, 
               lambda = 1., 
               optim = TRUE, 
               use_mono = FALSE,
               n_points = 10,
               maxit = 300, 
               eps_f = 1e-5, 
               eps_g = 1e-6,
               tol_simplex = 1e-10,
               precision = 1e-10,
               breaks = 100)
  names_ctrl <- names(ctrl)
  names_control <- names(control)
  ctrl[names_control] <- control
  if (length(unknown <- names_control[!names_control %in% names_ctrl])) {
    warning("unknown names in control: ", paste(unknown, collapse = ", "))
  }
  
  if (method == "stratified") {
    if (!is.null(delta)) {
      if (ctrl$eps_d == 0. || ctrl$eps_d > (1. - delta)/delta) {
        stop("control$eps_d should be strictly between 0 and (1. - delta)/delta")
      }
    }
    if (ctrl$use_mono & ctrl$optim) {
      warning("argument control$use_mono ignored when control$optim = TRUE")
    }
    set_ignored <- c("precision", "breaks")
  } else if (method == "vanilla") {
    set_ignored <- c("eps_d", "optim", "use_mono", "n_points", 
                     "maxit", "eps_f", "eps_g", "tol_simplex",
                     "precision", "breaks")
  } else {
    set_ignored <- c("eps_d", "optim", "use_mono", "n_points", 
                     "maxit", "eps_f", "eps_g", "tol_simplex")
  }
  
  if (length(ignored <- names_control[names_control %in% set_ignored])) {
    warning("argument ignored in control: ", paste(ignored, collapse = ", "))
  }
  out <- cpp_rnegmix(n, par, family, method, ctrl)
  out[!names(out) %in% c("pairs")] <- lapply(out[!names(out) %in% c("pairs")], as.numeric)
  return(out)
}
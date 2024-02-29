pnegmix <- function(q, par, family = c("normal", "gamma")) {
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
  
  if (family == "normal") {
    return(as.numeric(pnormix(q, par$w_p, par$mean_p, par$sd_p) - pnormix(q, par$w_n, par$mean_n, par$sd_n)))
  } else if (family == "gamma") {
    return(as.numeric(pgammix(q, par$w_p, par$shape_p, 1./par$rate_p) - pgammix(q, par$w_n, par$shape_n, 1./par$rate_n)))
  }
}
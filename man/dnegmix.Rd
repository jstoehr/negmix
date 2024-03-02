\name{dnegmix}
\alias{dnegmix}
\alias{rnegmix}
\alias{qnegmix}
\alias{pnegmix}
\title{Signed mixtures of distributions}
\description{
Probability density function, cumulative distribution function, quantile function and random generator for a mixture of distributions with signed weights.  
}
\usage{
dnegmix(x, par, family = c("normal", "gamma")) 
pnegmix(q, par, family = c("normal", "gamma"))
qnegmix(p, par, family = c("normal", "gamma"), 
        precision = 1e-10, breaks = 100)
rnegmix(n, par, family = c("normal", "gamma"), 
        delta = NULL, 
        method = c("stratified", "vanilla", "inv_cdf"), 
        control = list())
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of draws from the model.}
  \item{par}{list of model parameters. See ``Details''.}
  \item{family}{family of distributions used in the mixture.}
  \item{precision}{numerical precison to compute the approximate inverse distribution function. A piecewise linear function is used.}
  \item{breaks}{initial number of interpolating points to be used.}
  \item{delta}{acceptance rate target for the stratified sampling scheme.}
  \item{method}{the sampling method to be used. See ``Details''.}
  \item{control}{list of control parameters. See ``Details''.}
}
\value{
  \code{dnegmix} returns the probability density function computed at \code{x}.

  \code{pnegmix} returns the cumulative distribution function computed at \code{q}.

  \code{qnegmix} returns an approximation of the quantile function computed at \code{p}.

  \code{rnegmix} generates a random sample and returns a list with various diagnostic arguments:
  \describe{
  \item{sample}{Vector of the generated random sample.}
  \item{time}{The computing time (in seconds) required to get the sample.}
  \item{prob_accept}{The theoretical acceptance probability of the accept-reject algorithm. Only provided when \code{method} is \code{"stratified"} or \code{"vanilla"}.}
  \item{freq_accept}{The observed acceptance frequency of the accept-reject algorithm. Only provided when \code{method} is \code{"stratified"} or \code{"vanilla"}.}
\item{pairs}{Matrix; by column the pairs (two-component signed mixture) found by the simplex method. First row is the indices of the positive weight components and second row is the indices for the associated negative weight components. Only provided when \code{method} is \code{"stratified"}.}
\item{w_pair_p, w_pair_n}{Vectors of the non-normalised weights for, respectively, the positive and negative component of each pairs. Only provided when \code{method} is \code{"stratified"}.}
\item{residuals_p, r_p, residuals_n, r_n}{Vectors for the indices and the weights of the residuals of the positive weight components, and the indices and the weights of the residuals of the negative weight components, respectively. Only provided when \code{method} is \code{"stratified"}.}
  }
  
}
\details{
  A signed mixture writes as
  \deqn{
  \displaystyle\sum_{k = 1}^P w_{p,k} f_{k}(x) - \sum_{k = 1}^N w_{n,k} g_{k}(x),
  \quad\text{where}\quad (f_k, g_k)\in\code{family}.
  }
  The \code{par} argument is a list that specifies model parameters according to the \code{family} argument:
  \describe{
  \item{w_p, w_n}{Vectors of weights for the positive and negative components, respectively. Should only contains positive values.}
    \item{mean_p, sd_p, mean_n, sd_n}{Vectors for the mean and standard deviation parameters of the positive weight components, and the mean and standard deviation parameters of the negative weight components, respectively, for a mixture of normal distributions (\code{family = "normal"}).}
    \item{shape_p, rate_p, shape_n, rate_n}{Vectors for the shape and rate parameters of the positive weight components, and the shape and rate parameters of the negative weight components, respectively, for a mixture of gamma distributions (\code{family = "gamma"}).}
  }
  
  The \code{method} argument corresponds to the following sampling methods:
  \describe{
  \item{stratified}{Accept-reject algorithm using piecewise proposal distributions (fastest method in general).}
  \item{vanilla}{Accept-reject algorithm using the mixture of the positive weight components of the signed mixture as a proposal distribution.}
  \item{inv_cdf}{Inverse sampling with a piecewise linear approximation of the inverse distribution function.}
  }
  
  The \code{control} argument is a list that can supply any of the following inputs. It is mainly useful to compute the piecewise proposal distributions (histograms for each pairs of components) for the accept-reject method when \code{method = "stratified"}.
  \describe{
  \item{n_points}{The initial number of cells for the histogram. Only used by \code{method = "stratified"}. Default to 10.}
  \item{use_mono}{Logical; if TRUE the histogram computation uses the monoticity of the mixture on cells where it is known; if FALSE the histogram is computed as a difference of histograms. Only used by \code{method = "stratified"}. Default to FALSE.}
  \item{optim}{Logical; if TRUE the location of all local extrema of the mixture are computed (using the L-BFGS method) and added to the initial breaking points. The histogram is then computed using that the mixture distribution is a monotonic function on each cell. Only used by \code{method = "stratified"}. Default to TRUE.}
  \item{eps_d}{Largest possible error on the histogram to achieve the acceptance rate \code{delta} on a given pair of components. Ranges strictly betwenn 0 and (1 - \code{delta})/\code{delta}. Only used by \code{method = "stratified"}. Default to 0.9 * (1. - \code{delta})/\code{delta}.}
  \item{maxit}{The maximum number of iterations for the L-BFGS method. Only used by \code{method = "stratified"}. Default to 300.}
  \item{eps_f}{Stopping rule parameter for the L-BFGS method. The algorithm stops if
  the relative tolerance for two consecutive function values is lower than \code{eps_f}. Only used by \code{method = "stratified"}. Default to 1e-5.}
  \item{eps_g}{Stopping rule parameter for the L-BFGS method. The algorithm stops if
  the norm of the gradient at \code{x} is lower than \code{eps_g} * max(1, ||x||). Only used by \code{method = "stratified"}. Default to 1e-6.}
  \item{tol_simplex}{The floating point tolerance to be used in tests of equality in the simplex method. Only used by \code{method = "stratified"}. Default to 1e-10.}
  \item{precision}{numerical precison to compute the approximate inverse distribution function. A piecewise linear function is used. Only used by \code{method = "inv_cdf"}. Default to 1e-10.}
  \item{breaks}{initial number of interpolating points to be used. Only used by \code{method = "inv_cdf"}. Default to 100.}
  }
}
\references{
  Robert, C.P. and Stoehr J. (2024) Simulating signed mixtures. \emph{arXiv}:2401.16828.
}
\examples{
# --- Example of signed mixture of normal distributions
par <- render_negmix(5, 6, "normal")
# --- CDF example
x_min <- min(c(par$mean_p, par$mean_n) - 2 * c(par$sd_p, par$sd_n))
x_max <- max(c(par$mean_p, par$mean_n) + 2 * c(par$sd_p, par$sd_n))
x <- seq(x_min, x_max, .001)
y <- pnegmix(x, par, family = "normal")
plot(x, y, type = "l", main = "Signed mixture cdf")
# --- Inverse CDF example
plot(
  x, qnegmix(y, par, family = "normal"),
  main = "Comparison between theoretical and approximate quantiles",
  pch = 20
)
abline(0, 1, col = "firebrick", lty = 2)
# --- Random generators with stratified scheme
n <- 1e3
z_s <- rnegmix(n, par,
  family = "normal", method = "stratified",
  delta = 0.5 + 0.5 / sum(par$w_p)
)
hist(z_s$sample, freq = FALSE, breaks = floor(n / 20))
curve(
  dnegmix(x, par, family = "normal"),
  col = "royalblue", add = TRUE
)
# --- Comparison of the acceptance rates for the
# --- accept-reject based methods
z_v <- rnegmix(n, par, family = "normal", method = "vanilla")
z_s$accept_rate / z_v$accept_rate
# --- Relative efficiency of the stratified scheme
# --- compared to the approximate inverse cdf
z_q <- rnegmix(n, par, family = "normal", method = "inv_cdf")
z_q$time / z_s$time
}
\author{\packageAuthor{negmix}
Maintainer: \packageMaintainer{negmix}}


\name{render_negmix}
\alias{render_negmix}
\title{Signed mixtures of distributions}
\description{
\code{render_negmix} is used to generate random examples of mixture of distributions with signed weights.  
}
\usage{
render_negmix(k_p, k_n, family = c("normal", "gamma"), 
              by_pair = FALSE, control = list())
}
\arguments{
  \item{k_p, k_n}{number of positive and negative weight components, respectively.}
  \item{family}{family of distributions used in the mixture.}
  \item{by_pair}{logical; if TRUE the model generated is a mixture of two-component signed mixtures and positive weight components (hence exhibits a natural pairing structure); if FALSE the model generated is a mixture of signed mixtures, each with a single negative component (hence does not exhibit a natural pairing structure)}
  \item{control}{list of control parameters. See ``Details''.}
}
\value{
  \code{render_negmix} returns a list with weights and parameters for the distributions associated with the \code{family} argument:
  \describe{
    \item{w_p, w_n}{Vectors of weights for the positive and negative components, respectively.}
    \item{mean_p, sd_p, mean_n, sd_n}{Vectors for the mean and standard deviation parameters of the positive weight components, and the mean and standard deviation parameters of the negative weight components, respectively, for a mixture of normal distributions (\code{family = "normal"}).}
    \item{shape_p, rate_p, shape_n, rate_n}{Vectors for the shape and rate parameters of the positive weight components, and the shape and rate parameters of the negative weight components, respectively, for a mixture of gamma distributions (\code{family = "gamma"}).}
  }
  
}
\details{
The function aims at producing an example of signed mixture of the form
  \deqn{
  \displaystyle\sum_{k = 1}^P w_{p,k} f_{k}(x) - \sum_{k = 1}^N w_{n,k} g_{k}(x),
  \quad\text{where}\quad (f_k, g_k)\in\code{family}.
  }
The method first initializes the parameters of each negative weight component \eqn{g_j}by associating it to a positive weight component \eqn{f_j} to form a two-component signed mixture
  \deqn{
  w^{\rm pair}_{p,j} f_{j}(x) - w^{\rm pair}_{n,j} g_{j}(x).
  }

The \code{control} argument is a list that can supply any of the following inputs:
  \describe{
  \item{p_min, p_max}{Lower and upper bound for \eqn{1\big/\sum_{k = 1}^P w_{p,k}}. Control how close the negative part of the signed mixture is to the positive part. Default to 0 and 0.95.}
  \item{p_star}{The probability that, when initializing a negative component, the two-component signed mixture computed is associated with the largest possible weight \eqn{w^{\rm pair}_{n,j}} for the negative component. Default to 0.8.}
  \item{p_common}{The probability that, when initializing a negative component, \eqn{f_{j}} and \eqn{g_{j}} share common parameters with the associated positive component. Default to 0.2.}
  \item{p_rank}{The probability that, when initializing a negative component, \eqn{1 - w^{\rm pair}_{n,j}\big/w^{\rm pair}_{p,j}} is between \code{p_min} and \code{p_max}. Default to 1.}
  \item{k_p_init}{positive integer. The \code{k_p_init} first positive components are randomly set; if \code{k_p_init}\eqn{\leq}\code{k_n}, the \code{k_p_init} first positive components are used to initialized the negative components parameters; otherwise the initialization is achieved with a random set of components among the \code{k_p_init} first. Default to a random value, constrained to be lower than \code{min(k_n, k_p)} when \code{by_pair = TRUE}.}
  \item{k_p_pair}{The number of positive weight components to be used to build the two-component signed mixtures included in the model. The \code{k_p_pair - k_p_init} positive components are initialized from the \code{k_n} negative components. Ignored when \code{by_pair = FALSE}. Ranges between \code{k_p_init} and \code{k_p}. Default to a random value.}
    \item{n_pair}{The minimum number of two-component signed mixture to include in the model. Ignored when \code{by_pair = FALSE}. Ranges between \code{k_n + control$k_p_pair - control$k_p_init} and \code{k_n * (control$k_p_pair - control$k_p_init + 1)}. Default to a random value.}
    \item{n_single}{The number of positive residuals to include in the model, that is positive components \eqn{f_k} such that for any negative component \eqn{g_k}, there exists no weight to define a two-component signed mixture between \eqn{f_k} and \eqn{g_k}. Ignored when \code{by_pair = FALSE}. Ranges between 0 and \code{k_p - k_p_pair}. Default to a random value.}
  \item{maxit_0}{The maximum number of iterations for the Newton-Raphson method (used to find 0 a signed mixture and ensure positivity). Default to 300.}
  \item{eps_0}{Stopping rule parameter for the Newton-Raphson method. The algorithm stops if the function value is less than \code{eps_0}. Default to 1e-10.}
  \item{maxit}{The maximum number of iterations for the L-BFGS method (used to find local optimum and correct default of positivity). Default to 300.}
  \item{eps_f}{Stopping rule parameter for the L-BFGS method. The algorithm stops if
  the relative tolerance for two consecutive function values is lower than \code{eps_f}. Default to 1e-5.}
  \item{eps_g}{Stopping rule parameter for the L-BFGS method. The algorithm stops if
  the norm of the gradient at \code{x} is lower than \code{eps_g} * max(1, ||x||). Default to 1e-6.}
  }
}
\references{
  Robert, C.P. and Stoehr J. (2024) Simulating signed mixtures. \emph{arXiv}:2401.16828.
}
\examples{
# --- Example of signed mixture of gamma distributions
par <- render_negmix(5, 10,
  family = "gamma", by_pair = TRUE,
  control = list(p_max = 0.1, n_pair = 22)
)
1. / sum(par$w_p)
# --- Example of signed mixture of normal distributions
fam <- "normal"
k_p <- sample(1:50, 1)
k_n <- sample(1:50, 1)
par <- render_negmix(k_p, k_n,
  family = fam,
  by_pair = FALSE,
  control = list(p_min = 1e-2, p_max = 5e-1, p_star = 1.)
)
# --- PDF of the signed mixture
x_min <- min(c(par$mean_p, par$mean_n) - 2 * c(par$sd_p, par$sd_n))
x_max <- max(c(par$mean_p, par$mean_n) + 2 * c(par$sd_p, par$sd_n))
plot(x <- seq(x_min, x_max, .001),
  dnegmix(x, par, family = fam),
  type = "l", main = "Signed mixture pdf"
)
}
\author{\packageAuthor{negmix}
Maintainer: \packageMaintainer{negmix}}


\name{render_negmix}
\title{Signed mixtures of distributions}
\description{
\code{render_negmix} is used to generate random examples of mixture of distributions with signed weights.  
}
\usage{
render_negmix(k_p, k_n, family = c("normal", "gamma"), 
              par = list(), by_pair = FALSE, control = list())
}
\arguments{
  \item{k_p, k_n}{number of positive and negative weight components, respectively.}
  \item{family}{family of distributions used in the mixture.}
  \item{par}{list of parameters to initialized the positive weight components. See ``Details''.}
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
  The \code{par} argument is a list that specifies initial parameters for some postive weight components, according to the \code{family} argument. Each of the \code{k_n} negative weight components is initialized so that it forms a two-component signed mixture with one of these initial postive component.
  \describe{
    \item{mean_p, sd_p}{Vectors for the mean and standard deviation parameters of the initial positive weight components (\code{family = "normal"}).}
    \item{shape_p, rate_p}{Vectors for the shape and rate parameters of the initial positive weight components (\code{family = "gamma"}).}
  }
  The \code{control} argument is a list that can supply any of the following inputs:
  \describe{
  \item{p_min, p_max}{Lower and upper bound for the inverse of the sum of the positive component weights. Control how close the negative part of the signed mixture is to the positive part.}
  \item{p_star}{The probability that, when initializing a negative component, the two-component signed mixture computed is associated with the largest possible weight for the negative component.}
  \item{p_common}{The probability that, when initializing a negative component, it shares common parameters with the associated positive component.}
  \item{p_rank}{The probability that, when initializing a negative component, the inverse of the positive component weight is between \code{p_min} and \code{p_max}.}
  \item{k_init}{The number of positive weight components used to initialize the negative weight components. When the argument \code{par} is provided and \code{by_pair = TRUE}, the value is automatically set to the length of the initial parameters. Otherwise, default to a random value.}
  \item{k_pair}{The number of positive weight components used to build two-component signed mixtures included in the model. Ignored when \code{by_pair = FALSE}. Should be between \code{k_init} and \code{k_p}. Default to a random value between \code{k_init} and \code{k_p}.}
    \item{n_pair}{The number of two-component signed mixture to include in the model. Ignored when \code{by_pair = FALSE}. Should be between \code{k_n + k_pair - k_init} and \code{k_p * k_n}. Default to \code{k_n + k_pair - k_init}.}
    \item{k_single}{The number of positive weight components that cannot be associated to any negative weight components of the model. Ignored when \code{by_pair = FALSE}. Should be between 0 and \code{k_p - k_pair}. Default to a random value between 0 and \code{k_p - k_pair}.}
  \item{maxit_0}{The maximum number of iterations for the Newton-Raphson method. Default to 300.}
  \item{eps_0}{Stopping rule parameter for the Newton-Raphson method. The algorithm stops if the function value is less than \code{eps_0}. Default to 1e-10.}
  \item{maxit}{The maximum number of iterations for the L-BFGS method. Default to 300.}
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
# --- Example of signed mixture of normal distributions
render_negmix(5, 6, "normal")
render_negmix(5, 6, "normal", 
)
}
\author{\packageAuthor{negmix}
Maintainer: \packageMaintainer{negmix}}


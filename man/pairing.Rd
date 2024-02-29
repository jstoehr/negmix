\name{pairing}
\title{Signed mixtures of distributions}
\description{
Simplex algorithm to find a decomposition of a signed mixture into a mixture of two-component signed mixture and a residual signed mixture.
}
\usage{
pairing(delta, par, family = c("normal", "gamma"), tol_simplex = 1e-10)
}
\arguments{
  \item{delta}{parameter of the simplex objective function. See ``Details''.}
  \item{par}{list of model parameters. See ``Details''.}
  \item{family}{family of distributions used in the mixture.}
  \item{tol_simplex}{The floating point tolerance to be used in tests of equality.}
}
\value{
  \code{pairing} returns a list that supplies the decomposition:
  \describe{
\item{pairs}{Matrix; by column the pairs (two-components signed mixture) found by the simplex method. First row is the indices of the positive weight components \eqn{f_k} and second row is the indices for the associated negative weight components \eqn{g_k}.}
\item{w_pair_p, w_pair_n}{Vectors of the non-normalised weights \eqn{w^{\rm pair}_{p, k}} and  \eqn{w^{\rm pair}_{n, k}} for, respectively, the positive and negative component of each pairs.}
\item{residuals_p, r_p, residuals_n, r_n}{Vectors for the indices and the weights \eqn{r_i} of the residuals of the positive weight components, and the indices and the weights \eqn{s_j} of the residuals of the negative weight components, respectively.}
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
  \item{\code{w_p}, \code{w_n}}{Vectors of weights for the positive and negative components, respectively. Should only contains positive values.}
    \item{\code{mean_p}, \code{sd_p}, \code{mean_n}, \code{sd_n}}{Vectors for the mean and standard deviation parameters of the positive weight components, and the mean and standard deviation parameters of the negative weight components, respectively, for a mixture of normal distributions (\code{family = "normal"}).}
    \item{\code{shape_p}, \code{rate_p}, \code{shape_n}, \code{rate_n}}{Vectors for the shape and rate parameters of the positive weight components, and the shape and rate parameters of the negative weight components, respectively, for a mixture of gamma distributions (\code{family = "gamma"}).}
  }
  
  The simplex method aims at finding a decomposition of the form
  \deqn{
\displaystyle\sum_{k=1}^K \{w^{\rm pair}_{p,k} f_k(x) - w^{\rm pair}_{n, k}g_k(x)\}
+\sum_{i = 1}^P r_i f_i(x) - \sum_{j = 1}^N s_j g_j(x),
\quad\text{such that}\quad 1 - \frac{w^{\rm pair}_{n, k}}{w^{\rm pair}_{p, k}} \leq \code{delta}.
  }
The argument \code{delta} controls the acceptance probability of the accept-reject methods for sampling from each two-component signed mixture \eqn{w^{\rm pair}_{p,k} f_k(x) - w^{\rm pair}_{n, k}g_k(x)} involved in the decomposition.
}
\references{
  Robert, C.P. and Stoehr J. (2024) Simulating signed mixtures. \emph{arXiv}:2401.16828.
}
\examples{
# --- Example of signed mixture of normal distributions
par <- render_negmix(5, 6, "normal", by_pair = FALSE)
delta <- 0.5 + 0.5 / sum(par$w_p)
pairing(delta, par, "normal")
}
\author{\packageAuthor{negmix}
Maintainer: \packageMaintainer{negmix}}


\name{negmix-package}
\alias{negmix-package}
\alias{negmix}
\docType{package}
\title{
\packageTitle{negmix}
}
\description{
\packageDescription{negmix}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{negmix}
\packageIndices{negmix}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{negmix}
Maintainer: \packageMaintainer{negmix}
}
\references{
Robert, C.P. and Stoehr J. (2024). Simulating signed mixtures. \emph{arXiv}:2401.16828.
}
\keyword{ package }

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
# --- Finding a pairing associated with a signed mixture
delta <- 0.5 + 0.5 / sum(par$w_p)
pairing(delta, par, "normal")
}

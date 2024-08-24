rm(list = ls())
library(negmix)
library(rbenchmark)
library(ggplot2)
library(dplyr)
library(kableExtra)
par(mai = c(0.8, 0.5, 0.5, 0.1), family = "serif", cex.lab = 1.2, cex.axis = 0.9)

# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------
# ------------------------------------------------------

case <- "norm"
if (case == "norm") {
  family <- "normal"
} else if (case == "gam") {
  family <- "gamma"
}
dir <- paste0("example/negmix_", case, "/")
simplex_r <- FALSE
save_plot <- FALSE
# ------------- EXAMPLE WITH NATURAL PAIRING ------------
ex <- 1
inflated <- TRUE
xp_name <- paste0("negmix_", case, "_ex", ex, "_inf", inflated * 1.)

# ------------- RANDOM EXAMPLE ------------
param <- 2
# setting <- sample(1:4, 1)
# seed <- sample(1:50, 1)
# set_p <- c(0., 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 2e-1, 3e-1)
# i <- sample(2:length(set_p), 1)
# p_min <- set_p[i -1]
# p_max <- set_p[i]
setting <- 1
seed <- 4
p_min <- 5e-2
p_max <- 1e-1
xp_name <- paste0(
  "rand_neg_", case, "_set", setting, "_param", param,
  "_seed", seed, "_pmin", p_min, "_pmax", p_max
)
# --- Mountain
# filename <- "example/negmix_norm/rand_neg_norm_set3_param1_seed24_pmin1e-04_pmax0.001.RData"

# ------------
filename <- paste0(dir, xp_name, ".RData")
load(filename)

if (case == "norm") {
  x_1 <- min(c(par$mean_p, par$mean_n) - 2 * c(par$sd_p, par$sd_n))
  x_n <- max(c(par$mean_p, par$mean_n) + 2 * c(par$sd_p, par$sd_n))
} else if (case == "gam") {
  x_1 <- 0
  x_n <- qgamma(.95, max(par$shape_p), min(par$rate_p))
}

# --- Target density
(k_p <- length(par$w_p))
(k_n <- length(par$w_n))

x <- seq(x_1, x_n, .001)
plot(
  x, dnegmix(x, par, family),
  type = "l", xlab = "", ylab = ""
)

plot(
  x, pnegmix(x, par, family),
  type = "l", xlab = "", ylab = ""
)

# --- out_vanilla acceptance probability
cat("out_vanilla acceptance probability: ", 1 / sum(par$w_p), "\n")

# --- Building pairs
delta <- .4

# --- Set of possible pairs
map <- cpp_map_pairs(delta, par, family)
sum(map > 0)
# --- Pairing found with simplex method
out <- pairing(delta, par, family)

m_pairs <- matrix(0, k_p, k_n)
for (k in 1:ncol(out$pairs)) {
  m_pairs[out$pairs[1, k], out$pairs[2, k]] <- out$w_pair_p[k]/out$w_pair_n[k]
}

par(mfrow = c(1, 3), oma = c(0, 0, 0, 1))
image(map_pairs_m, axes = FALSE)
image(map, axes = FALSE)
image(m_pairs, axes = FALSE)
par(mfrow = c(1, 1))

n <- 1e5
# --- Global reject
out_vanilla <- rnegmix(n, par, family, method = "vanilla")

# --- Accept reject with pairs (with histogram)
out_strat <- rnegmix(n, par, family, delta)

# --- Numerical inverse of the cdf
out_inv <- rnegmix(n, par, family, method = "inv_cdf")

# --- Checking output sample
out <- out_strat
pl <- ggplot(data = data.frame(y = out$sample), aes(x = y))
pl <- pl + theme_bw() +
  theme(legend.position = c(.87, .87)) +
  labs(x = element_blank(), y = element_blank())
pl <- pl + geom_histogram(aes(y = after_stat(density)),
                          binwidth = (max(out$sample) - min(out$sample))/50,
                          fill = "grey70", colour = "grey70", alpha = .8)
pl <- pl + stat_function(fun = dnegmix,
                         args = list(par, family),
                         n = 1000, lwd = 0.9, aes(colour = "m"), alpha = .8)
pl <- pl + scale_colour_manual(name = 'Legend', values = "royalblue")
pl
if (save_plot) {
  ggsave(filename = paste0(dir, "figures/", "hist_", xp_name, ".jpeg"),
         width = 7, dpi = 300)
  # ggsave(filename = paste0(dir, "figures/", "hist_", xp_name, ".pdf"),
  #        width = 7)
}

# --- Comparison of acceptance probability
out_vanilla$prob_accept
out_vanilla$freq_accept

out_strat$prob_accept
out_strat$freq_accept

# --- Comparison of execution time
out_vanilla$time
out_strat$time
out_inv$time

# Relative efficiency
out_vanilla$time / out_strat$time
out_inv$time / out_strat$time
out_vanilla$time / out_inv$time

benchmark(
  rnegmix(n, par, family, method = "vanilla"),
  rnegmix(n, par, family, delta),
  rnegmix(n, par, family, method = "inv_cdf"),
  columns = c("test", "replications", "relative", "elapsed"),
  order = "elapsed", replications = 10
)

# --- Generate results for a single example
set.seed(2024)
n <- 1e2
delta <- c(.4, .6, .8)
set_eps <- c(.1, .2, .5, 1.)

out_vanilla <- rnegmix(n, par, family, method = "vanilla")
out_inv <- rnegmix(n, par, family, method = "inv_cdf", control = list(breaks = 500))

dat <- data.frame(rep(0, 4))
for (d in delta) {
  eps <- set_eps[set_eps < (1 - d)/d]
  for (e in eps) {
    out_strat <- rnegmix(n, par, family, d, control = list(eps_d = e))
    dat <- cbind(data.frame(dat,
                            res = c(
                              out_strat$prob_accept,
                              out_strat$freq_accept,
                              out_vanilla$time/out_strat$time,
                              out_inv$time/out_strat$time
                            )
    )
    )
  }
}

dat <- dat[, -1]
dat <- cbind(dat, data.frame(vanilla = c(out_vanilla$prob_accept, out_vanilla$freq_accept, 1,
                                         out_inv$time/out_vanilla$time),
                             row.names = c("prob", "freq", "R", "S")))

dat
dat %>%
  kbl(caption = "Example",
      format = "latex",
      row.names = TRUE,
      col.names = c("eps1", "eps2", "eps3", "eps4",
                    "eps1", "eps2", "eps3", "eps1", "eps2", "Vanilla"),
      digits = 3,
      align="r") %>%
  kable_minimal(full_width = F,  html_font = "Source Sans Pro")

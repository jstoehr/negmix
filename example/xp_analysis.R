rm(list = ls())
library(viridis)
library(ggplot2)
library(scales)
library(gridExtra)

save_fig <- FALSE
# --- All model 
dir_fig <- "example/figures/"
if (!dir.exists(dir_fig)) {
  dir.create(dir_fig, recursive = TRUE)
}
set_case <- c("norm", "gam")

# --- Single model
family <- "norm"
dir_fig_fam <- paste0("example/negmix_", family, "/figures/")
if (!dir.exists(dir_fig_fam)) {
  dir.create(dir_fig_fam, recursive = TRUE)
}
file_fam <- paste0(dir_fig_fam, "rand_negmix_", family)

# --- Loading result
size <- 500
res <- NULL
for (case in set_case) {
  source("xp_combine.R")
  
  if (case == "norm") {
    output$model <- "Normal signed mixture"
  } else if (case == "gam") {
    output$model <- "Gamma signed mixture"
  }
  
  if (is.null(res)) {
    res <- output
  } else {
    res <- rbind(res, output)
  }
}

n_res <- nrow(res)
res$cat[res$cat == 8] <- 7 # Changing the last category to 35%
res$r <- res$time_ar/res$time
res$q <- res$time_inv/res$time
res$q_ar <- res$time_inv/res$time_ar
res$prop_k <- res$k_p/(res$k_p + res$k_n) * 100
res$prop_pair <- res$n_pairs/res$n_max * 100

res$n <- as.factor(res$n)
res$cat <- as.factor(res$cat)
res$delta <- as.factor(res$delta)
res$eps <- as.factor(res$eps)
res$cat_pair <- as.factor(res$n_max %/% size)
res$model <- factor(res$model, levels = c("Normal signed mixture", 
                                          "Gamma signed mixture"))

temp <- c(0, (as.numeric(levels(res$cat_pair)) + 1) * size)
labels_pair <- rep(0, length(temp) - 1)
for (i in 1:(length(temp) - 1)) {
  labels_pair[i] <- paste0("(", temp[i], ", ", temp[i + 1],"]")
}

my_labels <- c("< 1e-4", "(1e-4, 0.001]", 
               "(0.001, 0.01]", "(0.01, 0.05]", 
               "(0.05, 0.10]", "(0.10, 0.20]", "(0.20, 0.35]")
my_labels <- c("< 0.01", "(0.01, 0.1]", 
               "(0.1, 1]", "(1, 5]", 
               "(5, 10]", "(10, 20]", "(20, 35]")

# --- Variety of examples (Focus on one model)
family <- "norm"
keep <- res$n == 10
output <- res[keep, ]
keep <- output$delta == .8
output <- output[keep, ]
keep <- output$eps == .1
output <- output[keep, ]
if (family == "norm") {
  keep <- output$model == "Normal signed mixture"
} else if (family == "gam") {
  keep <- output$model == "Gamma signed mixture"
}
output <- output[keep, ]
# --- Frequence for each subset of proba
data.frame(matrix(table(output$cat[output$model == "Normal signed mixture"]), 1))
data.frame(matrix(table(output$cat[output$model == "Gamma signed mixture"]), 1))

# --- Distribution of k_p and k_n
n_out <- nrow(output)

cat <- rep(output$cat, 2)
k <- c(output$k_p, output$k_n)
Components <- c(rep("Positive", n_out), rep("Negative", n_out))
data <- data.frame(cat, k, Components)

ggplot(data, aes(x = as.factor(cat), y = k, fill = Components)) + geom_boxplot() +
  theme_bw() +
  theme(legend.position="top",
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) + 
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam, "_components_number.pdf"), width = 7, height = 3)
  ggsave(filename = paste0(file_fam, "_components_number.jpeg"), width = 7, height = 3, dpi = 300)
}

# --- Distribution of the proportion
ggplot(output, aes(x = cat, y = prop_k)) + geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13)) +
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels)

if (save_fig) {
  # ggsave(filename = paste0(file_fam, "_proportion_components.pdf"), width = 7, height = 3)
  ggsave(filename = paste0(file_fam, "_proportion_components.jpeg"), width = 7, height = 3, dpi = 300)
}

# --- Distribution of the proportion of pairs selected
ggplot(output, aes(x = cat, y = n_max)) + geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13)) +
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels)

if (save_fig) {
  # ggsave(filename = paste0(file_fam, "_max_n_pairs.pdf"), width = 7, height = 3)
  ggsave(filename = paste0(file_fam, "_max_n_pairs.jpeg"), width = 7, height = 3, dpi = 300)
}

ggplot(output, aes(x = cat, y = prop_pair)) + theme_bw() +
  geom_boxplot() +
  labs(x = element_blank(), y = element_blank()) +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13)) +
  scale_x_discrete(labels = my_labels)

if (save_fig) {
  # ggsave(filename = paste0(file_fam, "_proportion_pairs_selected.pdf"), width = 7, height = 3)
  ggsave(filename = paste0(file_fam, "_proportion_pairs_selected.jpeg"), width = 7, height = 3, dpi = 300)
}

# --- Variety of examples (Focus on all models)
keep <- res$n == 10
output <- res[keep, ]
n_out <- nrow(output)

x <- rep(output$cat, 5)
# ---
y <- c(output$k_p, output$k_n, output$prop_k, output$n_max, output$prop_pair)
Legend <- c(rep("Positive", n_out), rep("Negative", n_out), rep(NA, 3 * n_out))
# ---
model <- rep(output$model, 5)
model <- factor(model, levels = c("Normal signed mixture", "Gamma signed mixture"))
# --- Plot title
lev_labels <- c("No. components",
                "% positive components", 
                "No. acceptable pairs", 
                "% pairs picked by simplex")
label <- factor(c(rep(lev_labels[1], 2 * n_out),
                  rep(lev_labels[2], n_out),
                  rep(lev_labels[3], n_out),
                  rep(lev_labels[4], n_out)), 
                levels = lev_labels)
# ---
data <- data.frame(x, y, Legend, model, label)

ggplot(data, aes(x = x, y = y, fill = Legend)) + geom_boxplot() +
  facet_grid(label ~ model, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 15), 
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 16)
  ) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  scale_fill_brewer(palette = "Oranges", na.translate = FALSE)

if (save_fig) {
  # ggsave(filename = paste0(dir_fig, "rand_negmix_variety.pdf"), width = 15, height = 12)
  ggsave(filename = paste0(dir_fig, "rand_negmix_variety.jpeg"), width = 15, height = 12, dpi = 300)
}

# --- Time efficiency (Focus on one model)
if (family == "norm") {
  keep <- res$model == "Normal signed mixture"
} else if (family == "gam") {
  keep <- res$model == "Gamma signed mixture"
}
output <- res[keep, ]
n_out <- nrow(output)

# --- Effect of eps and delta on the running time of the stratified scheme
ggplot(output, aes(x = delta, y = time, fill = eps)) + geom_boxplot() + 
  theme_bw() + 
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_y_continuous(trans = 'log10', labels = scientific) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam, "_impact_of_delta_eps.pdf"),
  #        width = 7, height = 3)
  ggsave(filename = paste0(file_fam, "_impact_of_delta_eps.jpeg"),
         width = 7, height = 3, dpi = 300)
}

ggplot(res, aes(x = delta, y = time, fill = eps)) + geom_boxplot() + 
  facet_wrap(~model, scales = "free_y") +
  theme_bw() + 
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 16),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_y_continuous(trans = 'log10', labels = scientific) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(dir_fig, "rand_negmix_impact_of_delta_eps.pdf"), 
  #        width = 15, height = 3)
  ggsave(filename = paste0(dir_fig, "rand_negmix_impact_of_delta_eps.jpeg"), 
         width = 15, height = 3, dpi = 300)
}

# --- Running time
data <- data.frame(
  cat = rep(res$cat, 2),
  time = c(res$time, res$time_ar),
  n = rep(res$n, 2),
  method = c(rep("Stratified method", n_out), 
             rep("Vanilla method", n_out))
)

ggplot(data, aes(x = cat, y = time, fill = n)) +  geom_boxplot() +
  facet_wrap(~method) +
  theme_bw() +
  theme(legend.position="bottom", 
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam, "_running_time.pdf"),
  #        width = 14, height = 3)
  ggsave(filename = paste0(file_fam, "_running_time.jpeg"),
         width = 14, height = 3, dpi = 300)
}

# --- Relative efficiency wrt to the vanilla acceptance probability
ggplot(output, aes(x = cat, y = r, fill = n)) + geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam, 
  #                          "_relative_efficiency_strat_compared_to_vanilla.pdf"),
  #        width = 7, height = 3)
  ggsave(filename = paste0(file_fam, 
                           "_relative_efficiency_strat_compared_to_vanilla.jpeg"),
         width = 7, height = 3, dpi = 300)
}


x <- c(res$cat, res$cat, res$cat)
y <- c(res$r, res$time, res$time_ar)
n <- c(res$n, res$n, res$n)
model <- c(res$model, res$model, res$model)
lev_labels <- c("Relative efficiency",
                "Stratified method", 
                "Vanilla method")
label <- factor(c(rep(lev_labels[1], n_res),
                  rep(lev_labels[2], n_res),
                  rep(lev_labels[3], n_res)), 
                levels = lev_labels)
data <- data.frame(x, y, n, model, label)

ggplot(data, aes(x = x, y = y, fill = n)) + geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  facet_grid(label~model, scales = "free_y") +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(dir_fig, "rand_negmix_time_efficiency.pdf"), 
  #        width = 15, height = 12)
  ggsave(filename = paste0(dir_fig, "rand_negmix_time_efficiency.jpeg"), 
         width = 15, height = 12, dpi = 300)
}


# --- Relative efficiency wrt to the choice of eps and delta
ggplot(output, aes(x = eps, y = r, fill = n)) + geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  facet_wrap(~delta, scales = "free_x") +
  theme_bw() +
  theme(legend.position="bottom", 
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam,
  #                          "_relative_efficiency_strat_compared_to_vanilla_per_delta_eps.pdf"),
  #        width = 7, height = 3)
  ggsave(filename = paste0(file_fam,
                           "_relative_efficiency_strat_compared_to_vanilla_per_delta_eps.jpeg"),
         width = 7, height = 3, dpi = 300)
}

ggplot(res, aes(x = eps, y = r, fill = n)) + geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  facet_grid(model~delta, scales = "free_x") +
  theme_bw() +
  theme(legend.position="bottom", 
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(dir_fig,
  #                          "rand_negmix_relative_efficiency_strat_compared_to_vanilla_per_delta_eps.pdf"),
  #        width = 15, height = 7)
  ggsave(filename = paste0(dir_fig,
                           "rand_negmix_relative_efficiency_strat_compared_to_vanilla_per_delta_eps.jpeg"),
         width = 15, height = 7, dpi = 300)
}


# --- Comparison of the relative efficiency of stratified and vanilla scheme compared to
# --- the numerical inverse of the cdf
data <- data.frame(
  n = rep(output$n, 2), 
  cat = rep(output$cat, 2), 
  q = c(output$q, output$q_ar), 
  method = c(rep("Stratified method", n_out), rep("Vanilla method", n_out))
)

ggplot(data, aes(x = cat, y = q, fill = n)) + geom_boxplot() + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  facet_wrap(~method) +
  theme_bw() +
  theme(legend.position="bottom", 
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam,
  #                          "_relative_efficiency_ar_method_compared_to_inv_cdf.pdf"),
  #        width = 14, height = 3)
  ggsave(filename = paste0(file_fam,
                           "_relative_efficiency_ar_method_compared_to_inv_cdf.jpeg"),
         width = 14, height = 3, dpi = 300)
}

x <- c(res$cat, res$cat)
y <- c(res$q, res$q_ar)
n <- c(res$n, res$n)
model <- c(res$model, res$model)
lev_labels <- c("Stratified method", 
                "Vanilla method")
label <- factor(c(rep(lev_labels[1], n_res),
                  rep(lev_labels[2], n_res)), 
                levels = lev_labels)
data <- data.frame(x, y, n, model, label)

ggplot(data, aes(x = x, y = y, fill = n)) + geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  facet_grid(label~model) +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(dir_fig, "rand_negmix_relative_efficiency_compared_to_inv_cdf.pdf"), 
  #        width = 15, height = 8)
  ggsave(filename = paste0(dir_fig, "rand_negmix_relative_efficiency_compared_to_inv_cdf.jpeg"), 
         width = 15, height = 8, dpi = 300)
}


# --- Comparison of the relative efficiency in fonction of the number of possible pairs
ggplot(output, aes(x = cat_pair, y = r, fill = n)) + geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  theme_bw() +
  theme(legend.position="bottom", 
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = labels_pair) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(file_fam,
  #                          "_relative_efficiency_compared_to_vanilla_per_pairs.pdf"),
  #        width = 14, height = 3)
  ggsave(filename = paste0(file_fam,
                           "_relative_efficiency_compared_to_vanilla_per_pairs.jpeg"),
         width = 14, height = 3, dpi = 300)
}


ggplot(res, aes(x = cat_pair, y = r, fill = n)) + geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  facet_wrap(~model, scales = "free_x") +
  theme_bw() +
  theme(legend.position="bottom", 
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = labels_pair) +
  scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(dir_fig, "rand_negmix_relative_efficiency_compared_to_vanilla_per_pairs.pdf"), 
  #        width = 15, height = 5)
  ggsave(filename = paste0(dir_fig, "rand_negmix_relative_efficiency_compared_to_vanilla_per_pairs.jpeg"), 
         width = 15, height = 5, dpi = 300)
}


ggplot(res, aes(x = cat, y = prob_s)) + geom_boxplot() +
  facet_grid(model~delta) +
  theme_bw() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 16)) + 
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(labels = my_labels) +
  #scale_y_continuous(trans = 'log10', labels = scientific, n.breaks = 6) +
  scale_fill_brewer(palette="Oranges")

if (save_fig) {
  # ggsave(filename = paste0(dir_fig, "rand_negmix_stratified_acceptance_rate.pdf"), 
  #        width = 15, height = 8)
  ggsave(filename = paste0(dir_fig, "rand_negmix_stratified_acceptance_rate.jpeg"), 
         width = 15, height = 8, dpi = 300)
}


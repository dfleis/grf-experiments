###################################################################################################################
# synthetic-split-stability.R
# Empirical demonstration of unstable gradient-based splits when the gradient of the local 
# expected score \nabla_{\theta,\nu} E[\psi_{\theta^*(x),\nu^*(x)}(O_i) | X_i = x] is nearly
# singular (and its estimator A_P).
#
#
#
# Generated 2025-Jan-05 23:59 EST
###################################################################################################################
library(ggplot2)
library(tidyr)
library(purrr)
library(dplyr)
library(grid)
library(gridExtra)
library(grf)

#-------------------------------------------------------
#----- HELPER FUNCTIONS
#-------------------------------------------------------
generate_data <- function(n, r, sigma = 1) {
  # n: Number of samples
  # r: Correlation between regressors W_1 and W_2
  X <- matrix(runif(n), nrow = n)
  W1 <- rnorm(n)
  W2 <- r * W1 + sqrt(1 - r^2) * rnorm(n)
  W <- cbind(W1, W2)
  
  eps <- rnorm(n, 0, sigma)
  Y <- (X[,1] > 0.5) * W[,1] + W[,2] + eps
  list(X = X, Y = Y, W = W)
}

lm_forest_stumps <- function(data, ...) {
  X <- data$X
  Y <- data$Y
  W <- data$W
  lm_forest(X = X, Y = Y, W = W, Y.hat = 0 * Y, W.hat = 0 * W, 
            min.node.size = nrow(X)/2 - 1,
            honesty = F, 
            alpha = 0, 
            ci.group.size = 1,
            compute.oob.predictions = F, ...)
}
stump_splits <- function(forest) {
  sapply(forest$`_split_values`, function(x) {
    if (length(x) > 3) stop("Tree is not a stump. Something went wrong.")
    if (length(x) < 3) stop("No split was made. Something went wrong.")
    return (x[1])
  })
}

#-------------------------------------------------------
#----- RUN SIMULATIONS
#-------------------------------------------------------
set.seed(1)
methods <- c("grad", "fpt1", "fpt2")
num.trees <- 2000
nreps <- 250
nobs <- 1000
rvals <- seq(0.8, 0.99, by = 0.01)# #seq(-0.98, 0.98, by = 0.02)

splits_list <- vector(mode = "list", length = length(rvals))
for (i in seq_along(rvals)) {
  r <- rvals[i]
  cat(sprintf("%s\tr = %0.4f\n", format(Sys.time()), r))
  
  sim <- lapply(1:nreps, function(rep_id) {
    data <- generate_data(nobs, r)
    
    res <- sapply(methods, function(method) {
      forest <- lm_forest_stumps(data, 
                                 method = method, 
                                 num.trees = num.trees,
                                 seed = rep_id) # same seed for every method
      stump_splits(forest)
    }, USE.NAMES = TRUE) # end sapply over method
    
    return (cbind(rep = rep_id, res))
  }) # end lapply over rep_id
  
  splits_list[[i]] <- do.call(rbind, sim) %>% 
    data.frame() %>%
    pivot_longer(cols = any_of(methods),
                 names_to = "Method",
                 values_to = "Split") %>%
    mutate(r = r)
}
# combine into single data frame
df <- do.call(rbind, splits_list) %>%
  mutate(Method = factor(Method, levels = methods)) %>%
  mutate(Method = recode(Method, fpt1 = "FPT1", fpt2 = "FPT2"))

# compute percentile statistics of the split variance
df_var_pctl <- df %>%
  group_by(Method, r, rep) %>% # compute split var per method x corr x rep
  summarise(var = var(Split)) %>%
  group_by(Method, r) %>% # compute percentiles aggregating over replications
  summarise(lo  = quantile(var, prob = 0.1),
            med = quantile(var, prob = 0.5),
            hi  = quantile(var, prob = 0.9))


#-------------------------------------------------------
#----- DRAW FIGURES
#-------------------------------------------------------
theme_common <- function(legend_position) {
  theme_minimal() + 
  theme(
    panel.border = element_rect(color = "black", fill = NA), 
    axis.ticks.length = unit(-0.2, "cm"),
    axis.ticks = element_line(),
    axis.line = element_line(),
    strip.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = legend_position,
    legend.title.align = 0.5,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(family = "mono"))
}


# variance median/percentile bands: splits vs. regressor correlation
# grad, FPT1 (FPT) -- one panel
plt_var <- df_var_pctl %>%
  filter(Method != "FPT2") %>%
  mutate(Method = recode(Method, FPT1 = "FPT")) %>%
  ggplot(aes(x = r, y = med, color = Method, fill = Method)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, linewidth = NA) + 
  geom_line(linewidth = 1) + 
  geom_point(size = 2) + 
  labs(x = "Regressor Correlation", y = "Split Variance", 
       title = sprintf("Median split variance over %s replications", nreps)) + 
  scale_y_log10(labels = function(x) sprintf("1e-%d", -log10(x))) + 
  theme_common(legend_position = c(0.9, 0.6))


# boxplots: splits vs. regressor correlation
# grad, FPT1 (FPT) -- vertically facetted
plt_box <- df %>%
  filter(Method != "FPT2") %>%
  filter(abs(r) >= 0.8) %>%
  mutate(Method = recode(Method, FPT1 = "FPT")) %>%
  ggplot(aes(x = as.factor(r), y = Split, color = Method, fill = Method)) + 
  geom_hline(yintercept = 0.5, color = "gray50", linewidth = 0.5, linetype = "dashed") +
  geom_boxplot(linewidth = 0.5, outlier.alpha = 0, alpha = 0.25, width = 0.5) + 
  facet_grid(rows = vars(Method)) + 
  labs(x = "Regressor Correlation", y = "Split Value", 
       title = "Boxplots: Split values over regressor correlations") + 
  scale_y_continuous(limits = c(0, 1)) +
  theme_common(legend_position = c(0.1, 0.5)) + 
  scale_x_discrete(breaks = pretty(rvals, 7))

# combine into single figure
g <- arrangeGrob(
  plt_box,
  plt_var,
  heights = unit(c(3.2, 2.8), "inch"),
  nrow = 2,
  padding = unit(0, "line"))

#-------------------------------------------------------
#----- SAVE FIGURE
#-------------------------------------------------------
filename <- function(ext) sprintf("figures/grf-synthetic-split-stability.%s", ext)
ggsave(filename("pdf"), g, width = 7.5, height = 6, units = "in")
ggsave(filename("png"), g, width = 7.5, height = 6, units = "in")


         
  
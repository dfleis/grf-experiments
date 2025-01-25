library(tidyverse)

source("utils-example.R")

### Setup/parameters
seed_data <- seed_grf <- 1

args_grf <- list(num.trees = 2000,
                 min.node.size = 5,
                 ci.group.size = 2,
                 seed = seed_grf)

treatment_funs <- list(function(x) 0 * x[,1], # baseline
                       function(x) 5 * (x[,1] - 0.2) * (x[,1] - 0.7) - 0.2,
                       function(x) 2 * (x[,1] > 0.5) - 1,
                       function(x) cos(4 * pi * x[,1]))

args_example <- list(model_type = "hte", 
                     n_test = 250, 
                     n_train = 2500, 
                     p = 1,
                     sigma_eps = 0.5,
                     effect_funs = treatment_funs,
                     args_grf = args_grf,
                     seed = seed_data)

### Fit model
df <- do.call(run_example, args_example)

### Make plots
format_arg <- function(arg_name, arg_list, default_val = NULL, arg_lab = NULL) {
  val <- ifelse(arg_name %in% names(arg_list), arg_list[[arg_name]], default_val)
  if (!is.null(arg_lab)) arg_name <- arg_lab
  paste(val, arg_name)
}

model_type_str <- ifelse(args_example$model_type == "hte", 
                         "heterogeneous treatment effects", 
                         "varying coefficient effects") 
title_str <- sprintf("Example: Estimating %s", model_type_str)
subtitle_str <- paste0("Data: ", 
                       args_example$n_test, " test covariates (sequential), ", 
                       args_example$n_train, " training samples\n", 
                       "Forest: ", 
                       format_arg("num.trees", args_grf, 2000, "trees"), ", ",
                       100 * attr(df, "ci_level"), " % CI")

# xx <- df$est$X1
# uxx <- sort(unique(xx))
# df_points <- df$est %>% filter(X1 %in% uxx[seq(1, length(uxx), length.out = 10)])

df$est %>%
  ggplot(aes(x = X1, y = value, color = method, fill = method)) + 
  labs(x = expression("Covariate"~X[1]), y = "Effect Value", color = "Method", fill = "Method",
       title = title_str, subtitle = subtitle_str) + 
  facet_grid(rows = vars(effect), cols = vars(method),
             labeller = label_parsed) + 
  geom_hline(yintercept = 0, color = "gray90", linewidth = 0.75) + 
  geom_line(data = df$truth, color = "gray30", linewidth = 0.4, linetype = "solid") + 
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), color = NA, alpha = 0.3) +
  geom_line(linewidth = 0.75, alpha = 0.9) +
  # geom_errorbar(data = df_points, aes(ymin = ci_lo, ymax = ci_hi), color = "black", alpha = 0.25) + 
  # geom_point(data = df_points, size = 0.75, color = "black") +
  theme_minimal() + 
  theme(
    panel.border = element_rect(fill = NA), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-0.2, "cm"),
    axis.ticks = element_line(),
    axis.line = element_line(),
    legend.position = "none",
    strip.text.x = element_text(family = "mono", size = 12))



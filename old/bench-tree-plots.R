suppressMessages(library(tidyverse))
library(gridExtra)
library(scales)
library(pals)

path <- "data/bench-tree"
#basename_pattern <- "bench-lm-stumps.*nrep10-niter50.csv"
basename_pattern <- "bench-lm-trees.*nrep10-niter50.csv"
#basename_pattern <- "bench-multi_arm-stumps.*nrep10-niter50.csv"


filenames <- list.files(path = path, pattern = basename_pattern, full.names = T)
df_all <- lapply(filenames, function(fnm) read.csv(fnm)) %>% 
  bind_rows() %>%
  mutate(method = recode(method, "grad" = "grad", "fpt1" = "FPT1", "fpt2" = "FPT2")) %>%
  mutate(method = factor(method, levels = c("grad", "FPT1", "FPT2")),
         K = as.factor(K),
         p = as.factor(p),
         n = as.factor(n))

# df_all <- df_all %>%
#   filter(p %in% c(1, 5))


#----------------------------------------------------------------------
#---------- RESTRUCTURE DATA FOR PLOTTING
#----------------------------------------------------------------------
# median aggregate over (non-GC) iterations
df_reps <- df_all %>%
  #filter(!gc) %>% # filter iterations that had GC
  group_by(method, K, p, n, rep) %>%
  summarise(time   = median(time),
            splits = median(splits), .groups = "keep")

# median aggregate over replications
df_summary <- df_reps %>%
  group_by(method, K, p, n) %>%
  summarise(time   = median(time),
            splits = mean(splits), .groups = "keep")

# relative performance measures: new method/old method
df_ratio <- df_summary %>%
  group_by(K, p, n) %>% {
    time_frame <- reframe(.,
                          "FPT1" = time[method == "FPT1"]/time[method == "grad"],
                          "FPT2" = time[method == "FPT2"]/time[method == "grad"]) %>%
      gather(method, time, "FPT1":"FPT2") 
    splits_frame <- reframe(., 
                            "FPT1" = splits[method == "FPT1"]/splits[method == "grad"],
                            "FPT2" = splits[method == "FPT2"]/splits[method == "grad"]) %>%
      gather(method, splits, "FPT1":"FPT2")
    by <- join_by(K, p, n, method)
    inner_join(time_frame, splits_frame, by = by)
  } %>% mutate(hline = 1)

#----------------------------------------------------------------------
#---------- PLOT CONFIG / HELPERS
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 12
MY_FONT_SIZE_LEGEND <- 12
MY_FONT_SIZE_LABEL <- 11
MY_SHAPES <- c(21, 22, 24)
MY_COLORS <- scales::hue_pal()(3)
MY_FILLS <- scales::hue_pal()(3)

my_theme <- function(legend_pos = NULL) {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.major.x = element_line(color = "gray95"),
    legend.position = "inside",
    legend.position.inside = legend_pos,
    legend.key = element_rect(color = NA),
    legend.key.width = unit(1.5, "lines"),
    legend.margin = margin(0.15, 0.5, 0.15, 0.3, "lines"),
    legend.background = element_rect(linewidth = 0.5, color = "gray65"),
    legend.title = element_blank(),
    #legend.title = element_text(hjust = 1),
    legend.text = element_text(size = MY_FONT_SIZE_LEGEND, 
                               family = MY_FONT_FAMILY_MONO),
    text = element_text(size = MY_FONT_SIZE, family = MY_FONT_FAMILY)) 
}
my_common_elements <- function(xlab = "Regressors (K)", ylab, ...) {
  list(xlab(xlab),
       ylab(ylab),
       labs(...),
       geom_line(linewidth = 0.5),
       geom_point(size = 2.5, stroke = 0.75),
       facet_grid(vars(n), vars(p)))
}
my_summary_elements <- function(...) {
  c(my_common_elements(...),
    list(scale_shape_manual(values = MY_SHAPES),
         scale_color_manual(values = MY_COLORS),
         scale_fill_manual(values = MY_FILLS)))
}
my_ratio_elements <- function(...) {
  list(geom_hline(aes(yintercept = hline), 
                  col = 'gray75', linewidth = 0.75, linetype = "solid"),
       my_common_elements(...),
       scale_shape_manual(values = MY_SHAPES[-1]),
       scale_color_manual(values = MY_COLORS[-1]),
       scale_fill_manual(values = MY_FILLS[-1]))
}
my_scale_y <- function(expand = c(0.06, 0.05), ...) {
  scale_y_continuous(expand = expand, ...)
}
my_scale_y_percent <- function(limits = c(0, 1.2),
                               expand = c(0.01, 0.01),
                               breaks = seq(0, 1, 0.25), ...) {
  my_scale_y(limits = limits,
             expand = expand,
             breaks = breaks,
             labels = scales::percent_format(accuracy = 1), ...)
}
my_labeller <- function(df, x, y) {
  if (missing(x)) x <- (1 + length(levels(df$K)))/2
  
  df_label <- df %>%
    select(method, n, p) %>%
    distinct() %>%
    group_by(method, n, p) %>%
    summarise(label   = paste0("n = ", n, ", p = ", p),
              n_label = paste0("n = ", n),
              p_label = paste0("p = ", p), .groups = "keep") %>%
    mutate(x = x, y = y)
  
  geom_label(data = df_label, 
             mapping = aes(x = x, y = y, label = label),
             size = MY_FONT_SIZE_LABEL, 
             size.unit = "pt",
             color = "gray30", 
             label.size = 0,
             fill = "white", 
             family = MY_FONT_FAMILY,
             label.r = unit(0, "lines"))
}

#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
# df_ratio %>%
#   ggplot(aes(x = K, y = time, shape = method, fill = method, color = method, group = method)) + 
#   my_ratio_elements(ylab = "Median runtime reduction") +
#   my_scale_y_percent() +
#   my_labeller(df_ratio, y = 1.15) + 
#   my_theme(legend_pos = c(0.15, 0.2))
df_ratio %>%
  ggplot(aes(x = K, y = 1/time, shape = method, fill = method, color = method, group = method)) +
  my_ratio_elements(ylab = "Median speedup factor", subtitle = "Varying coefficient model") +
  ggtitle("Timing benchmarks (stumps): Fit time speedup factor over GRF-grad") +
  my_scale_y(expand = c(0.1, 0.1)) + 
  my_labeller(df_ratio, y = 0.25) +
  my_theme(legend_pos = c(0.1, 0.93))

# df_summary %>%
#   ggplot(aes(x = K, y = time*1000, shape = method, fill = method, color = method, group = method)) +
#   my_summary_elements(ylab = "Runtime (ms)") +
#   scale_y_log10() +
#   my_labeller(df_ratio, y = max(df_summary$time*1000)) +
#   my_theme(legend_pos = c(0.08, 0.88))

# df_reps %>%
#   ggplot(aes(x = K, y = time, color = method, fill = method)) +
#   xlab("Regressor features (K)") +
#   ylab("Runtime (ms)") +
#   geom_violin(position = position_identity()) +
#   facet_grid(vars(p), vars(n)) +
#   my_theme(legend_pos = c(0.2, 0.85)) +
#   #geom_point(data = df_summary, mapping = aes(x = K, y = time)) +
#   geom_line(data = df_summary, mapping = aes(x = K, y = time, group = method)) +
#   scale_fill_manual(values = ggplot2::alpha(MY_FILLS, 0.1)) +
#   scale_color_manual(values = ggplot2::alpha(MY_FILLS, 0.8)) +
#   my_labeller(df_ratio, y = max(df_reps$time))
#   scale_y_continuous(trans = log2_trans(),
#                      breaks = trans_breaks("log2", function(x) 2^x),
#                      labels = trans_format("log2", math_format(2^.x)))
# 
# df_reps %>%
#   ggplot(aes(x = K, y = time, color = method, fill = method)) +
#   xlab("Regressor features (K)") +
#   ylab("Runtime (ms)") +
#   geom_point() +
#   facet_grid(vars(p), vars(n)) +
#   my_theme(legend_pos = c(0.2, 0.85))


# dtmp <- df_reps %>%
#   filter(K == 128)
# ggplot(dtmp, aes(x = method, y = time, color = method)) + 
#   geom_boxplot()
# 



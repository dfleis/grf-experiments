suppressMessages(library(tidyverse))
library(ggh4x)
library(scales)
library(pals)

path <- "data/OLD/bench-forest-3-3/"
basename_pattern_time <- "bench-forest-3-3-TIME-vcm-.*.csv"
basename_pattern_splits <- "bench-forest-3-3-SPLITS-vcm-.*.csv"

load_data <- function(basename_pattern) {
  filenames <- list.files(path = path, pattern = basename_pattern, full.names = T)
  lapply(filenames, function(fnm) read.csv(fnm)) %>% 
    bind_rows() %>%
    #filter(!gc) %>% # filter iterations that had GC
    mutate(method = recode(method, "grad" = "grad", "fpt1" = "FPT1", "fpt2" = "FPT2")) %>%
    mutate(method = factor(method, levels = c("grad", "FPT1", "FPT2")),
           K = as.factor(K),
           p = as.factor(p),
           n = as.factor(n),
           num.trees = as.factor(num.trees),
           setting_id = as.factor(setting_id))
}

df_init <- list(time = basename_pattern_time, 
               splits = basename_pattern_splits) %>%
  lapply(load_data)

group_vars <- c("method", "model_type", "setting_id", "num.trees", "K", "p", "n")
measure_vars <- c("time", "splits")
gc_vars <- c("gc", "gc0", "gc1", "gc2")

#-------------------------------------------------------
#----- Restructure data for ggplot
#-------------------------------------------------------
# Drop garbage collection indicator columns
df_all <- df_init %>% lapply(select, -any_of(gc_vars))

### Aggregate over bench iterations (for times) and trees (for split counts)
# Time: Median aggregate over bench iterations (repeat timing done over the same data)
#       but keep the different replications (repeat timing done over new samples of the data).
# Splits: Since iterations use the same data/seed, forests are identical over bench
#         iterations (given a fixed replication, data parameters, and GRF parameters).
#         For this reason, we don't track the iteration-level split data.
#         However, the raw split data contains the split count for each tree in the
#         forest. If we median aggregate the splits over all trees then we can pair
#         the median times (over iterations) with the median tree split count.
# Iterations are only used for timing to account for unexpected/difficult to control
# system noise that may affect the benchmark timing.
my_aggregate_fun <- function(x, extra_group_vars = NULL, FUN = median) {
  x %>%
    group_by(across(any_of(c(group_vars, extra_group_vars)))) %>% 
    summarise(across(any_of(measure_vars), FUN), .groups = "keep")
}
df_reps <- df_all %>% 
  modify_at("time", my_aggregate_fun, "rep") %>%   # aggregate times over iterations
  modify_at("splits", my_aggregate_fun, "rep") %>% # aggregate split count over trees
  reduce(inner_join, by = c(group_vars, "rep"))    # merge into a single data frame

### Aggregate over bench replications
# Time: Median aggregate over bench replications.
# Splits: Median aggregate over bench replications.
df_summary <- df_reps %>% my_aggregate_fun()

### Relative performance measures (grad/FPT)
my_ratio_fun <- function(x) {
  group_vars_ratio <- setdiff(group_vars, "method")
  
  calc_method_ratio <- function(target) { # helper
    function(x) {
      # pick(method) returns a tibble with the single column "method"
      m <- pick(method)[[1]] # extract the column vector
      x[m == "grad"] / x[m == target]
    }
  }
  
  # Calculate the GRF-grad/GRF-FPT ratio of the measure_vars variables/columns
  # over the different group_vars_ratio combinations (parameters/settings values).
  df_ratio_wide <- x %>%
    group_by(across(any_of(group_vars_ratio))) %>%
    summarise(across(any_of(measure_vars),
                     list("FPT1" = calc_method_ratio("FPT1"),
                          "FPT2" = calc_method_ratio("FPT2"))),
              .groups = "drop")
  
  # The original ratio data frame is wide-formatted with separate columns for the
  # GRF-grad/GRF-FPT1 and GRF-grad/GRF-FPT2 ratios. To use with ggplot, we want to
  # make the data long-formatted with one column for the ratio value and another
  # column for the ratio method/variable, i.e. a column whose entries indicate
  # whether the corresponding ratio value is grad/FPT1 or grad/FPT2.
  df_ratio_long <- df_ratio_wide %>%
    pivot_longer(cols = contains(measure_vars),
                 names_to = c("metric", "method"),
                 names_pattern = "(.*)_(.*)",
                 values_to = "value")
  
  # If we're calculating the ratio over several measure_vars at the same time,
  # then we want to separate the different ratio measurements into their own
  # columns, i.e. making the data frame slightly less long-formatted.
  df_ratio_long %>%
    pivot_wider(names_from = metric,
                values_from = value) %>%
    mutate(hline = 1) # add hline for use with ggplot (draw the baseline ratio)
}

# Actually calculate the ratios 
df_ratio <- df_summary %>% my_ratio_fun()

### Create data frames that drop FPT1 and only include grad and FPT2
drop_and_recode <- function(x) {
  x %>%
    filter(method %in% c("grad", "FPT2")) %>%
    mutate(method = recode(method, "FPT2" = "FPT"))
}

df_reps2 <- df_reps %>% drop_and_recode()
df_summary2 <- df_summary %>% drop_and_recode()
df_ratio2 <- df_ratio %>% drop_and_recode()

#----------------------------------------------------------------------
#---------- Plot setup/config/helpers
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 12
MY_FONT_SIZE_LEGEND <- 10
MY_FONT_SIZE_LABEL <- 10
MY_SHAPES <- c(21, 22, 24)
MY_COLORS <- scales::hue_pal()(3)
MY_FILLS <- scales::hue_pal()(3)

my_theme <- function(...) {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.major.x = element_line(color = "gray95"),
    text = element_text(size = MY_FONT_SIZE, family = MY_FONT_FAMILY),
    ...) 
}
my_theme_legend_common <- function(...) {
  theme(
    legend.key = element_rect(color = NA),
    legend.key.width = unit(1.5, "lines"),
    legend.margin = margin(0.15, 0.5, 0.15, 0.3, "lines"),
    legend.background = element_rect(linewidth = 0.5, color = "gray65"),
    legend.text = element_text(size = MY_FONT_SIZE_LEGEND,
                               family = MY_FONT_FAMILY), ...)
}
my_theme_legend_inside <- function(position = NULL, ...) {
  my_theme_legend_common(
    legend.position = "inside",
    legend.position.inside = position,
    ...
  )
}


my_labeller <- function(df, x, y) {
  if (missing(x)) x <- (1 + length(levels(df$K)))/2
  
  df_label <- df %>%
    select(method, n, p, num.trees, setting_id) %>%
    distinct() %>%
    group_by(method, n, p, num.trees, setting_id) %>%
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
simple_labeller <- function(labels) label_both(labels, sep = " = ")
combined_np_labeller <- as_labeller(function(x) {
  sapply(x, function(val) {
    parts <- strsplit(as.character(val), "\\.")[[1]]
    paste0("(n, p): (", parts[1], ", ", parts[2], ")")
  })
})

#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
MODEL_TYPE <- "vcm"
SETTING_ID <- 1
# lapply(df_all, function(df) unique(df$num.trees))
NUM_TREES <- 100
subtitle_str <- sprintf("%s Setting %s", toupper(MODEL_TYPE), SETTING_ID)
subtitle_str_trees <- sprintf("%s, num.trees = %s", subtitle_str, NUM_TREES)

df_ratio %>%
  filter(num.trees == NUM_TREES, model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = K, y = time, shape = method, color = method, group = method)) +
  labs(x = "Regressors (K)", y = "Speedup factor", shape = "Method", color = "Method") + 
  ggtitle("Fit time speedup factor: GRF-grad/GRF-FPT", subtitle = subtitle_str_trees) +
  geom_hline(aes(yintercept = hline), 
             col = 'gray75', linewidth = 0.75, linetype = "solid") + 
  geom_line(linewidth = 0.5) + 
  geom_point(size = 2.5, stroke = 0.75) + 
  facet_grid(vars(n), vars(p)) +
  scale_y_continuous(expand = c(0.1, 0.1)) + 
  my_labeller(df_ratio, y = 0.60) +
  my_theme()

df_ratio2 %>%
  filter(model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = K, y = time, color = num.trees, group = num.trees)) +
  labs(x = "Regressors (K)", y = "Speedup factor", shape = "num.trees", color = "num.trees") + 
  ggtitle("Fit time speedup factor: GRF-grad/GRF-FPT", subtitle = subtitle_str) +
  geom_hline(aes(yintercept = hline), 
             col = 'gray75', linewidth = 0.75, linetype = "solid") + 
  geom_line(linewidth = 0.5) + 
  geom_point(size = 2, stroke = 0.75, shape = 22) + 
  facet_grid(vars(n), vars(p)) +
  scale_y_continuous(expand = c(0.1, 0.1)) + 
  my_labeller(df_ratio, y = 0.60) +
  my_theme()

df_ratio2 %>%
  filter(num.trees == NUM_TREES, model_type == MODEL_TYPE) %>%
  ggplot(aes(x = K, y = time, color = setting_id, shape = setting_id, group = setting_id)) +
  labs(x = "Regressors (K)", y = "Speedup factor", shape = "Setting", color = "Setting") + 
  ggtitle("Fit time speedup factor: GRF-grad/GRF-FPT", subtitle = sprintf("num.trees = %s", NUM_TREES)) + 
  geom_hline(aes(yintercept = hline), 
             col = 'gray75', linewidth = 0.75, linetype = "solid") + 
  geom_line(linewidth = 0.5) + 
  geom_point(size = 2.5, stroke = 0.75) + 
  facet_grid(vars(n), vars(p)) +
  scale_y_continuous(expand = c(0.1, 0.1)) +
  my_labeller(df_ratio, y = 0.60) +
  my_theme()

df_summary %>%
  filter(num.trees == NUM_TREES, model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  ggtitle("Absolute fit times", subtitle = subtitle_str_trees) +
  labs(x = "", y = "Median fit time (sec)") + 
  geom_col(color = "black") +
  facet_nested_wrap(vars("(n, p)" = interaction(n, p), K), 
                    ncol = nlevels(df_summary$K),
                    scales = "free", 
                    nest_line = element_line(color = "gray75"),
                    labeller = labeller("(n, p)" = combined_np_labeller, "K" = label_both)) +
  theme(strip.background = element_blank(),
        axis.ticks = element_line("gray75"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "gray75", fill = NA),
        legend.position = "none")

df_summary2 %>%
  filter(num.trees == NUM_TREES, model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  ggtitle("Absolute fit times", subtitle = subtitle_str_trees) +
  labs(x = "", y = "Median fit time (sec)") + 
  geom_col(color = "black") +
  facet_nested_wrap(vars("(n, p)" = interaction(n, p), K), 
                    ncol = nlevels(df_summary2$K) * nlevels(df_summary2$n),
                    scales = "free", 
                    nest_line = element_line(color = "gray75"),
                    labeller = labeller("(n, p)" = combined_np_labeller, "K" = label_both)) +
  theme(strip.background = element_blank(),
        axis.ticks = element_line("gray75"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "gray75", fill = NA),
        legend.position = "none")


df_ratio2 %>%
  filter(model_type == MODEL_TYPE) %>%
  ggplot(aes(x = K, y = time, color = setting_id, shape = setting_id, group = setting_id)) +
  ggtitle("Fit time speedup factor: GRF-grad/GRF-FPT", subtitle = toupper(MODEL_TYPE)) +
  labs(color = "Setting", shape = "Setting") + 
  geom_hline(aes(yintercept = hline), 
             col = 'gray75', linewidth = 0.75, linetype = "solid") + 
  geom_point(size = 2.5, stroke = 0.5) + 
  geom_line() +
  facet_nested_wrap(vars(num.trees, "(n, p)" = interaction(n, p)), 
                    nrow = nlevels(df_ratio$num.trees),
                    #scales = "free", 
                    nest_line = element_line(color = "gray75"),
                    labeller = labeller("(n, p)" = combined_np_labeller, "num.trees" = label_both)) +
  scale_y_continuous(limits = range(0.8, df_ratio$time), expand = c(0.1, 0.1)) + 
  theme(strip.background = element_blank(),
        axis.ticks = element_line("gray75"),
        panel.grid.major.x = element_line(color = "gray95", linewidth = 0.5),
        panel.grid.major.y = element_line(color = "gray95", linewidth = 0.5),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "gray75", fill = NA))


df_summary2 %>%
  filter(model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  ggtitle("Absolute fit times", subtitle = subtitle_str) +
  labs(x = "", y = "Median fit time (sec)") + 
  geom_col(color = "black") +
  facet_nested_wrap(vars(num.trees, "(n, p)" = interaction(n, p), K), 
                    #ncol = nlevels(df_summary$time$K),
                    ncol = nlevels(df_summary$K) * nlevels(df_summary$n),
                    scales = "free", 
                    nest_line = element_line(color = "gray75"),
                    labeller = labeller("(n, p)" = combined_np_labeller, 
                                        "num.trees" = label_both,
                                        "K" = label_both)) + 
  theme(strip.background = element_blank(),
        axis.ticks = element_line("gray75"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(color = "gray75", fill = NA),
        legend.position = "none")

# 
# df_summary %>%
#   group_by(setting_id, num.trees, K, p, n, model_type, method) %>%
#   summarize(mean(splits)) %>%
#   filter(setting_id == 1, p == 2, num.trees == 1, K == 4)

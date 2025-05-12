suppressMessages(library(tidyverse))
library(ggh4x)
library(scales)
library(pals)

path <- "data/OLD/bench-tree-3-3"
basename_pattern <- "bench-tree-3-3-vcm-.*.csv"

filenames <- list.files(path = path, pattern = basename_pattern, full.names = T)

df_init <- lapply(filenames, function(fnm) read.csv(fnm)) %>% 
  bind_rows() %>%
  #filter(!gc) %>% # filter iterations that had GC
  mutate(method = recode(method, "grad" = "grad", "fpt1" = "FPT1", "fpt2" = "FPT2")) %>%
  mutate(method = factor(method, levels = c("grad", "FPT1", "FPT2")),
         K = as.factor(K),
         p = as.factor(p),
         n = as.factor(n),
         setting_id = as.factor(setting_id))

group_vars <- c("method", "model_type", "setting_id", "stump", "K", "p", "n")
measure_vars <- c("time", "splits")
gc_vars <- c("gc", "gc0", "gc1", "gc2")

#-------------------------------------------------------
#----- Restructure data for ggplot
#-------------------------------------------------------
# Drop garbage collection indicator columns
df_all <- df_init %>% select(-any_of(gc_vars))

### Aggregate over bench iterations
# Median aggregate over bench iterations (repeat timing done over the same data)
# but keep the different replications (repeat timing done over new samples of the data).
df_reps <- df_all %>%
  group_by(across(all_of(c(group_vars, "rep")))) %>% 
  summarise(across(all_of(measure_vars), median), .groups = "keep")

### Aggregate over bench replications
# Median aggregate over bench replications.
df_summary <- df_reps %>%
  group_by(across(all_of(group_vars))) %>% 
  summarise(across(all_of(measure_vars), median), .groups = "keep")

### Relative performance measures (grad/FPT)
calc_method_ratio <- function(target) {
  function(x) {
    # pick(method) returns a tibble with the single column "method"
    m <- pick(method)[[1]] # extract the column vector
    x[m == "grad"] / x[m == target]
  }
}
df_ratio <- df_summary %>%
  group_by(across(all_of(setdiff(group_vars, "method")))) %>%
  summarise(across(all_of(measure_vars), 
                   list("FPT1" = calc_method_ratio("FPT1"),
                        "FPT2" = calc_method_ratio("FPT2"))),
            .groups = "drop") %>%
  pivot_longer(cols = contains(measure_vars),
               names_to = c("metric", "method"),
               names_pattern = "(.*)_(.*)",
               values_to = "value") %>%
  pivot_wider(names_from = metric,
              values_from = value) %>%
  mutate(hline = 1) # Add hline to more easily plot the baseline ratio

#----------------------------------------------------------------------
#---------- Plot setup/config/helpers
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
SETTING_ID <- 3
STUMP <- FALSE

stump_str <- ifelse(STUMP, "Stumps (single splits)", "Trees")
subtitle_str <- sprintf("%s: %s Setting %s", stump_str, toupper(MODEL_TYPE), SETTING_ID)

df_ratio %>%
  filter(stump == STUMP, model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = K, y = time, shape = method, fill = method, color = method, group = method)) +
  my_ratio_elements(ylab = "Speedup factor", subtitle = subtitle_str) +
  ggtitle("Fit time speedup factor: GRF-grad/GRF-FPT") +
  my_scale_y(expand = c(0.1, 0.1)) + 
  my_labeller(df_ratio, y = 0.25) +
  my_theme(legend_pos = c(0.1, 0.93))

df_summary %>%
  filter(stump == STUMP, model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = method, y = time, fill = method)) +
  ggtitle("Absolute fit times", subtitle = subtitle_str) +
  labs(x = "", y = "Median fit time (sec)") + 
  geom_col(color = "black") +
  facet_nested_wrap(vars("(n, p)" = interaction(n, p), K), 
                    ncol = nlevels(df_all$K),
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

df_reps %>%
  filter(stump == STUMP, model_type == MODEL_TYPE, setting_id == SETTING_ID) %>%
  ggplot(aes(x = method, y = splits, fill = method)) +
  ggtitle("Split count", subtitle = subtitle_str) +
  labs(x = "", y = "Splits") +
  geom_boxplot(outliers = FALSE) +
  #geom_violin() +
  facet_nested_wrap(vars("(n, p)" = interaction(n, p), K), 
                    ncol = nlevels(df_all$K),
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


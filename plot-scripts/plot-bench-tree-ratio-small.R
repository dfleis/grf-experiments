library(tidyverse)
library(rlang)
library(ggh4x)
library(scales)
library(pals)

path <- "data/bench-tree-small-1-10/"
basename_pattern <- "bench-tree-small-1-10-.*.csv"

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
#---------- PLOT CONFIG / HELPERS
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 7
MY_FONT_SIZE_STRIP_X <- 8
MY_FONT_SIZE_STRIP_Y <- 8
MY_FONT_SIZE_LEGEND <- 8
MY_FONT_SIZE_LABEL <- 7
MY_FONT_SIZE_AXIS_X <- 9
MY_FONT_SIZE_AXIS_Y <- 8
MY_COLORS <- scales::hue_pal()(5) # number of settings

my_theme <- function() {
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = MY_FONT_SIZE_STRIP_X, family = MY_FONT_FAMILY_MONO, face = "bold"),
    strip.text.y = element_text(size = MY_FONT_SIZE_STRIP_Y),
    axis.ticks.x = element_line(color = "gray75", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(-0.05, "cm"),
    axis.text.x = element_text(family = MY_FONT_FAMILY, size = MY_FONT_SIZE_AXIS_X),
    axis.text.y = element_text(family = MY_FONT_FAMILY, size = MY_FONT_SIZE_AXIS_Y),
    axis.line.y = element_line(color = "gray90", linewidth = 0.5),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(color = "gray75", fill = NA),
    #panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.25),
    text = element_text(size = MY_FONT_SIZE, family = MY_FONT_FAMILY),
    legend.key = element_rect(color = NA),
    legend.key.width = unit(1, "lines"),
    legend.margin = margin(0.15, 0.5, 0.15, 0.3, "lines"),
    #legend.background = element_rect(linewidth = 0.5, color = "gray65"),
    legend.text = element_text(size = MY_FONT_SIZE_LEGEND))
}

custom_labeller <- function(name = "", sep = " = ") {
  function(x) paste0(name, sep, x)
}

my_labeller_p <- function(df, x, y) {
  if (missing(x)) x <- (1 + length(levels(df$K)))/2
  
  df_label <- df %>%
    select(p, K, setting_id) %>%
    distinct() %>%
    group_by(p, K, setting_id) %>%
    summarise(label   = paste0("dim(X) = ", p),
              .groups = "keep") %>%
    mutate(x = x, y = y)
  
  geom_label(data = df_label, 
             mapping = aes(x = x, y = y, label = label),
             size = MY_FONT_SIZE_LABEL, 
             size.unit = "pt",
             color = "gray50", 
             label.size = 0,
             fill = "white", 
             family = MY_FONT_FAMILY,
             label.r = unit(0, "lines"))
}

#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
MODEL_TYPE <- "hte"

STUMP <- c(TRUE, FALSE)
K_FILTER <- c(4, 16) # K = 4, 16
n_FILTER <- c("1000", "4000") # n = 1000, 2000, 4000

df_plt <- df_ratio %>%
  filter(
    model_type == MODEL_TYPE,
    stump %in% STUMP,
    K %in% K_FILTER,
    n %in% n_FILTER
  ) %>%
  mutate(stump2 = as.character(stump),
         stump2 = recode(stump2, "TRUE" = "Single split (stump)", "FALSE" = "Single tree")) %>%
  mutate(stump2 = factor(stump2, levels = c("Single split (stump)", "Single tree"))) %>%
  mutate(
    n = droplevels(n), 
    K = droplevels(K), 
    p = droplevels(p),
    setting_id = droplevels(setting_id),
    stump2 = droplevels(stump2)
  )

title_str <- "Fit time speedup factor: GRF-grad/GRF-FPT (single split/single tree)"

#----- MODEL-SPECIFIC SETTINGS
my_ylim <- switch(MODEL_TYPE,
                  "vcm" = c(0.62, 2.25),
                  "hte" = c(0.62, 2.25))
my_breaks <- pretty(seq(min(my_ylim), max(my_ylim), length.out = 5))

my_label_p_ypos <- switch(MODEL_TYPE, "vcm" = 2.1, "hte" = 2.1)
my_model_colors <- switch(MODEL_TYPE, "vcm" = MY_COLORS[1:4], "hte" = MY_COLORS)

subtitle_str <- 
  switch(MODEL_TYPE,
         "vcm" = "Varying coefficient model (VCM)",
         "hte" = "Heterogeneous treatment effects (HTE)")
#subtitle_str <- paste0("Single split & single tree\n", subtitle_str)
legend_str <- sprintf("%s\nSetting", toupper(MODEL_TYPE))

#----- MAKE PLOTS (barplot)
df_plt_adjusted <- df_plt %>%
  mutate(ratio_from_one = time - 1) # calculate height from baseline of 1


plt_bar <- df_plt_adjusted %>%
  ggplot(aes(x = K, y = ratio_from_one, fill = setting_id)) + 
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, col = 'gray75', linewidth = 0.75, linetype = "solid") +
  labs(x = "Regressor dimension (K)", y = "Speedup factor", fill = legend_str) +
  ggtitle(title_str, subtitle = subtitle_str) + 
  ggh4x::facet_nested(
    stump2 + n ~ method,
    axes = "x", remove_labels = "x",
    labeller = labeller(
      n = custom_labeller(name = "n", sep = " = ")),
    nest_line = element_line(color = "gray75")) +
  scale_y_continuous(
    oob = rescale_none,
    limits = my_ylim - 1,
    labels = function(x) sprintf("%.1f", x + 1),
    breaks = my_breaks - 1) +
  scale_fill_manual(values = my_model_colors) +
  my_labeller_p(df_plt, y = my_label_p_ypos - 1) +
  my_theme()
plt_bar

filename_plt_bar <- sprintf("figures/tree/bench-tree-ratio-small-%s.pdf", MODEL_TYPE)
ggsave(filename_plt_bar, plot = plt_bar, width = 4.5, height = 4.5)
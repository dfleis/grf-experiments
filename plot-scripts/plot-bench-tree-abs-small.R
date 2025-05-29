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

#----------------------------------------------------------------------
#---------- PLOT CONFIG / HELPERS
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 9
MY_FONT_SIZE_STRIP <- 9
MY_FONT_SIZE_LEGEND <- 12
MY_FONT_SIZE_AXIS_X <- 8
MY_FONT_SIZE_AXIS_Y <- 7
MY_COLORS <- MY_FILLS <- c("#f8766d", "#4390e8", "#4cd461") # reddish, blueish, greenish

my_theme <- function(legend_position = "none") {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = MY_FONT_SIZE_STRIP),
    axis.ticks.x = element_line(color = "gray75", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(-0.05, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    axis.text.x = element_text(size = MY_FONT_SIZE_AXIS_X, family = MY_FONT_FAMILY_MONO, face = "bold"),
    axis.text.y = element_text(size = MY_FONT_SIZE_AXIS_Y),
    axis.line.y = element_line(color = "gray90", linewidth = 0.5),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(color = "gray75", fill = NA),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.25),
    text = element_text(size = MY_FONT_SIZE, family = MY_FONT_FAMILY),
    legend.position = legend_position,
    plot.title = element_text(size = MY_FONT_SIZE, 
                              margin = margin(b = 0)),
    plot.subtitle = element_text(size = MY_FONT_SIZE, 
                                 margin = margin(b = -4)))
}

#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
MODEL_TYPE <- "vcm"
SETTING_ID <- 4

STUMP <- c(TRUE, FALSE)
K_FILTER <- c(4, 16) # K = 4, 16
n_FILTER <- c("1000", "4000") # n = 1000, 2000, 4000

df_plt <- df_summary %>%
  filter(
    model_type == MODEL_TYPE, 
    setting_id == SETTING_ID,
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
    setting_id = droplevels(setting_id)
  )


#----- MODEL-SPECIFIC SETTINGS
model_label <- switch(MODEL_TYPE, 
                      "vcm" = "VCM", 
                      "hte" = "HTE")
title_str <- sprintf("Median fit times: %s (single split/single tree)", model_label)
subtitle_str <- sprintf("%s Setting %s", toupper(MODEL_TYPE), SETTING_ID)

#----- MAKE PLOTS (barplot)
plt <- df_plt %>%
  ggplot(aes(x = method, y = time, fill = method)) + 
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Method", y = "Fit time (seconds)", fill = "Method") +
  ggtitle(title_str, subtitle = subtitle_str) +
  ggh4x::facet_nested(
    stump2 + n ~ K, 
    scales = "free",independent = "y", 
    #nrow = 2,
    labeller = labeller(
      #nt = as_labeller(function(x) paste("Trees =", x)),
      n = as_labeller(function(x) paste("n =", x)),
      K = as_labeller(function(x) paste("K =", x))),
    nest_line = element_line(color = "gray75"),
    strip = strip_nested(size = "variable")) + 
  scale_fill_manual(values = MY_COLORS) + 
  my_theme()
plt

filename_plt <- sprintf("figures/tree/bench-tree-abs-small-%s-%s.pdf", MODEL_TYPE, SETTING_ID)
ggsave(filename_plt, plot = plt, width = 4, height = 4.5)
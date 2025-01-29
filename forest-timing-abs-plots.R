suppressMessages(library(tidyverse))
library(ggh4x)
library(gridExtra)
library(scales)
library(pals)

path <- "data/forest-timing"
basename_pattern_hte <- "forest-timing-hte-.*-long-.*\\.csv"
basename_pattern_vcm <- "forest-timing-vcm-.*-long-.*\\.csv"

filenames_hte <- list.files(path = path, pattern = basename_pattern_hte, full.names = T)
filenames_vcm <- list.files(path = path, pattern = basename_pattern_vcm, full.names = T)

df_list_hte <- lapply(filenames_hte, function(fnm) read.csv(fnm))
df_list_vcm <- lapply(filenames_vcm, function(fnm) read.csv(fnm))

df_all <- list(df_list_hte, df_list_vcm) %>% # absolute time
  bind_rows() %>%
  mutate(method = recode(method, "grad" = "grad", "fpt2" = "FPT")) %>%
  mutate(method = factor(method, levels = c("grad", "FPT")),
         K = as.factor(K),
         p = as.factor(p),
         n = as.factor(n),
         nt = as.factor(nt),
         setting_id = as.factor(setting_id),
         model_type = as.factor(model_type))

#----------------------------------------------------------------------
#---------- PLOT CONFIG / HELPERS
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 8
MY_FONT_SIZE_AXIS_X <- 7
MY_FONT_SIZE_AXIS_Y <- 5
MY_FONT_SIZE_STRIP <- 6
MY_COLORS <- scales::hue_pal()(2) # methods

my_theme <- function() {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = MY_FONT_SIZE_STRIP,
                              margin = margin(2,0,2,0)),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    axis.text.x = element_text(size = MY_FONT_SIZE_AXIS_X),
    axis.text.y = element_text(size = MY_FONT_SIZE_AXIS_Y),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.major.x = element_line(color = "gray95"),
    text = element_text(size = MY_FONT_SIZE, family = MY_FONT_FAMILY),
    legend.position = "none",
    plot.title = element_text(size = MY_FONT_SIZE,
                              margin = margin(b = 0)),
    plot.subtitle = element_text(size = MY_FONT_SIZE,
                                 margin = margin(b = -4)))
}

#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
MODEL_TYPE <- "vcm"
SETTING_ID <- 1

df_plt <- df_all %>% filter(model_type == MODEL_TYPE, setting_id == SETTING_ID)

### MODEL-SPECIFIC SETTINGS
my_ylim <- switch(MODEL_TYPE,
                  "vcm" = c(0.6, 3.15),
                  "hte" = c(0.85, 1.65))
my_breaks <- pretty(seq(min(my_ylim), max(my_ylim), length.out = 5))

my_label_np_ypos <- switch(MODEL_TYPE, "vcm" = 0.75, "hte" = 0.90)
my_label_trees_ypos <- switch(MODEL_TYPE, "vcm" = 3.05, "hte" = 1.60)
my_model_colors <- switch(MODEL_TYPE, "vcm" = MY_COLORS[1:4], "hte" = MY_COLORS)

model_label <- switch(MODEL_TYPE,
                      "vcm" = "Varying coefficient model",
                      "hte" = "Heterogeneous treatment effects")
title_str <- sprintf("Forest fit times: %s", model_label)
subtitle_str <- sprintf("Setting %s", SETTING_ID)

### MAKE PLOTS
plt <- df_plt %>%
  ggplot(aes(x = method, y = median, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Method", y = "Fit time (seconds)", fill = "Method") +
  ggtitle(title_str, subtitle = subtitle_str) +
  facet_nested_wrap(
    nt ~ n + K,
    scales = "free_y",
    nrow = 3,
    labeller = labeller(
      nt = as_labeller(function(x) paste("Trees =", x)),
      n = as_labeller(function(x) paste("n =", x)),
      K = as_labeller(function(x) paste("K =", x))),
    nest_line = element_line(linewidth = 0.1, color = "gray50"),
    strip = strip_nested(size = "variable")) +
  my_theme()

filename_plt <- sprintf("figures/forest-timing-abs-%s-%s.png", MODEL_TYPE, SETTING_ID)
ggsave(filename_plt, plot = plt, width = 7, height = 4)

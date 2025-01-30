suppressMessages(library(tidyverse))
library(ggh4x)
library(gridExtra)
library(scales)
library(pals)

path <- "data/forest-timing-large-K"
basename_pattern_hte <- "forest-timing-large-K-hte-.*-long-.*\\.csv"
basename_pattern_vcm <- "forest-timing-large-K-vcm-.*-long-.*\\.csv"

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

MY_FONT_SIZE <- 10
MY_FONT_SIZE_AXIS_X <- 9
MY_FONT_SIZE_AXIS_Y <- 7
MY_FONT_SIZE_STRIP <- 8
#MY_COLORS <- scales::hue_pal()(2) # methods
MY_COLORS <- c("#4390e8", "#4cd461")

my_theme <- function() {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = MY_FONT_SIZE_STRIP, 
                              margin = margin(2,0,2,0)),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "gray90", linewidth = 0.25),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    axis.text.x = element_text(size = MY_FONT_SIZE_AXIS_X),
    axis.text.y = element_text(size = MY_FONT_SIZE_AXIS_Y),
    axis.line.y = element_line(color = "gray90", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.25),
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
MODEL_TYPE <- "hte"

df_plt <- df_all %>% filter(model_type == MODEL_TYPE)

### MODEL-SPECIFIC SETTINGS
model_label <- switch(MODEL_TYPE, 
                      "vcm" = "Varying coefficient model (VCM)", 
                      "hte" = "Heterogeneous treatment effects (HTE)")
title_str <- sprintf("Forest fit times: %s", model_label)
subtitle_str <- "K = 256, n = 10000"

### MAKE PLOTS
plt <- df_plt %>%
  ggplot(aes(x = method, y = median, fill = method)) + 
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Method", y = "Fit time (seconds)", fill = "Method") +
  ggtitle(title_str, subtitle = subtitle_str) +
  facet_nested_wrap(
    nt ~ setting_id, 
    scales = "free", 
    nrow = 4,
    labeller = labeller(
      nt = as_labeller(function(x) paste("Trees =", x)),
      n = as_labeller(function(x) paste("n =", x)),
      K = as_labeller(function(x) paste("K =", x)),
      setting_id = as_labeller(function(x) sprintf("%s Setting %s", toupper(MODEL_TYPE), x))),
    nest_line = element_line(linewidth = 0.1, color = "gray50"),
    strip = strip_nested(size = "variable")) + 
  scale_fill_manual(values = MY_COLORS) + 
  my_theme()

filename_plt <- sprintf("figures/forest-timing-large-K-abs-%s.png", MODEL_TYPE)
ggsave(filename_plt, plot = plt, width = 7, height = 5)

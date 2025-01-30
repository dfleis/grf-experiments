suppressMessages(library(tidyverse))
library(ggh4x)
library(gridExtra)
library(scales)
library(pals)

path <- "data/forest-mse/K4"
basename_pattern_hte <- "forest-mse-avg_mse-hte-.*\\.csv"
basename_pattern_vcm <- "forest-mse-avg_mse-vcm-.*\\.csv"

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

MY_FONT_SIZE <- 12
MY_FONT_SIZE_AXIS_X <- 10
MY_FONT_SIZE_AXIS_Y <- 8
MY_FONT_SIZE_STRIP <- 10
MY_COLORS <- MY_FILLS <- c("#4390e8", "#4cd461")

my_theme <- function() {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = MY_FONT_SIZE_STRIP, 
                              margin = margin(2,0,2,0)),
    axis.ticks.x = element_line(color = "gray90", linewidth = 0.25),
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
    legend.position = "none")
    #legend.key = element_rect(color = NA),
    #legend.key.width = unit(1, "lines"),
    #legend.margin = margin(0.15, 0.5, 0.15, 0.3, "lines"),
    #legend.background = element_rect(linewidth = 0.5, color = "gray65"),
    #legend.text = element_text(size = MY_FONT_SIZE_LEGEND))
}


#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
MODEL_TYPE <- "hte"
K_VAL <- 4
N_VAL <- 20000

df_plt <- df_all %>% filter(model_type == MODEL_TYPE, K == K_VAL, n == N_VAL)

### MODEL-SPECIFIC SETTINGS
model_label <- switch(MODEL_TYPE, 
                      "vcm" = "Varying coefficient model (VCM)", 
                      "hte" = "Heterogeneous treatment effects (HTE)")
title_str <- sprintf("MSE estimates: %s", model_label)
subtitle_str <- sprintf("5000 test observations, 40 replications\nK = %s, n = %s", K_VAL, N_VAL)


plt <- df_plt %>%
  ggplot(aes(x = method, y = avg_mse, fill = method, color = method)) + 
  labs(x = "Method", y = "MSE", fill = "Method", color = "Method") +
  ggtitle(title_str, subtitle = subtitle_str) +
  geom_boxplot(alpha = 0.3) + 
  #geom_boxplot(alpha = 0.3, outliers = FALSE, width = 0.5) + 
  #geom_violin(alpha = 0.3) +
  #geom_point(position = position_jitter(seed = 1, width = 0.2), 
  #           alpha = 0.5, size = 0.75) +
  facet_nested_wrap(
    . ~ setting_id,  
    scales = "free", 
    labeller = labeller(
      setting_id = as_labeller(function(x) sprintf("%s Setting %s", toupper(MODEL_TYPE), x)))) +
  scale_fill_manual(values = MY_COLORS) +
  scale_color_manual(values = MY_COLORS) +
  my_theme()

filename_plt <- sprintf("figures/forest-mse-%s-%s-%s.png", MODEL_TYPE, K_VAL, N_VAL)
ggsave(filename_plt, plot = plt, width = 7, height = 4)



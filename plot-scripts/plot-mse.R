library(tidyverse)
library(rlang)
library(ggh4x)
library(scales)
library(pals)

path <- "data/forest-mse/"
basename_pattern <- "forest-mse-avg_mse-.*.csv"

filenames <- list.files(path = path, pattern = basename_pattern, full.names = T)

df_init <- lapply(filenames, function(fnm) read.csv(fnm)) %>% 
  bind_rows() %>%
  #filter(!gc) %>% # filter iterations that had GC
  mutate(method = recode(method, "grad" = "grad", "fpt1" = "FPT1", "fpt2" = "FPT2")) %>%
  mutate(method = factor(method, levels = c("grad", "FPT1", "FPT2")),
         K = as.factor(K),
         p = as.factor(p),
         n = as.factor(n),
         num.trees = as.factor(nt),
         setting_id = as.factor(setting_id)) %>%
  select(-nt, -FUN, -effect_type, -prob_type)

#----------------------------------------------------------------------
#---------- PLOT CONFIG / HELPERS
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 9
MY_FONT_SIZE_STRIP <- 8
MY_FONT_SIZE_LEGEND <- 12
MY_FONT_SIZE_AXIS_X <- 8
MY_FONT_SIZE_AXIS_Y <- 7
MY_COLORS <- MY_FILLS <- c("#f8766d", "#4390e8", "#4cd461") # reddish, blueish, greenish

my_theme <- function(legend_position = "bottom") {
  legend_position <- match.arg(legend_position, choices = c("none", "left", "right", "top", "bottom", "inside"))
  thm <- theme(
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
    legend.position = legend_position
  ) 
  if (legend_position == "none") {
    return (thm)
  } else {
    thm <- thm +
      theme(
        legend.key = element_rect(color = NA),
        legend.key.width = unit(1, "lines"),
        legend.margin = margin(0.15, 0.5, 0.15, 0.3, "lines"),
        #legend.background = element_rect(linewidth = 0.5, color = "gray65"),
        legend.text = element_text(
          size = MY_FONT_SIZE_LEGEND, 
          family = MY_FONT_FAMILY_MONO,
          face = "bold"
        )
      )
    return (thm)
  }
}


custom_labeller <- function(name = "", sep = " = ") {
  function(x) paste0(name, sep, x)
}
null_labeller <- function() function(x) ""

#----------------------------------------------------------------------
#---------- DRAW PLOTS
#----------------------------------------------------------------------
MODEL_TYPE <- "vcm"
K_FILTER <- 4 # K = 4, 256

n_FILTER <- c("10000", "20000", "100000") # n = 10000, 20000, 100000
num.trees_FILTER <- c("100") # num.trees = 100

df_plt <- df_init %>%
  filter(
    model_type == MODEL_TYPE,
    K %in% K_FILTER,
    n %in% n_FILTER,
    num.trees %in% num.trees_FILTER
  ) %>%
  mutate(
    n = droplevels(n), 
    K = droplevels(K), 
    p = droplevels(p), 
    setting_id = droplevels(setting_id), 
    num.trees = droplevels(num.trees)
  )

#----- MODEL-SPECIFIC SETTINGS
model_label <- switch(
  MODEL_TYPE, 
  "vcm" = "Varying coefficient model (VCM)", 
  "hte" = "Heterogeneous treatment effects (HTE)"
)
title_str <- sprintf("MSE estimates: %s", model_label)
subtitle_str <- paste0("50 model replications, 5000 test observations\nK = ", K_FILTER)

plt_mse <- df_plt %>%
  mutate(avg_mse = 100 * avg_mse) %>%
  ggplot(aes(x = method, y = avg_mse, fill = method)) + 
  labs(x = "Method", y = expression(100%*%"MSE"), fill = "Method") + 
  ggtitle(title_str, subtitle = subtitle_str) + 
  geom_boxplot(alpha = 0.65, linewidth = 0.3, staplewidth = 0.25, outliers = F) +
  #geom_boxplot(alpha = 0.65, linewidth = 0.3, staplewidth = 0.25, outlier.alpha = 0.2) + 
  ggh4x::facet_nested(
    vars(num.trees, n), vars(setting_id),  
    scales = "free",
    independent = "y",
    axes = "x", remove_labels = "x",
    labeller = labeller(
      setting_id = custom_labeller(sprintf("%s Setting", toupper(MODEL_TYPE)), " "),
      n = custom_labeller("n", " = "),
      num.trees = custom_labeller("nTrees", " = ")
    ),
    nest_line = element_line(color = "gray75")
  ) +
  scale_fill_manual(values = MY_COLORS) +
  scale_color_manual(values = MY_COLORS) +
  my_theme("none") 
plt_mse


filename_plt <- sprintf("figures/forest/mse-%s-K%s.pdf", MODEL_TYPE, K_FILTER)
ggsave(filename_plt, plot = plt_mse, width = 7, height = 3.5)








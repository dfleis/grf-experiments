suppressMessages(library(tidyverse))
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
         
df_ratio <- df_all %>% # ratio of GRF-grad/GRF-FPT (speedup factor)
  group_by(model_type, setting_id, n, p, K, nt) %>%
  reframe("FPT" = median[method == "grad"]/median[method == "FPT"]) %>%
  gather(method, median, "FPT") %>%
  rename(ratio = median) %>%
  mutate(hline = 1)

#----------------------------------------------------------------------
#---------- PLOT CONFIG / HELPERS
#----------------------------------------------------------------------
# MY_FONT_FAMILY <- "Helvetica-Narrow"
# MY_FONT_FAMILY_MONO <- "monospace"
MY_FONT_FAMILY <- "sans"
MY_FONT_FAMILY_MONO <- "mono"

MY_FONT_SIZE <- 10
MY_FONT_SIZE_LEGEND <- 9
MY_FONT_SIZE_LABEL <- 10
MY_COLORS <- scales::hue_pal()(5) # number of settings

my_theme <- function() {
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.ticks.x = element_line(color = "gray90", linewidth = 0.25),
    axis.ticks.y = element_line(color = "gray90", linewidth = 0.25),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    axis.line.y = element_line(color = "gray90", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.25),
    panel.grid.major.x = element_line(color = "gray95", linewidth = 0.25),
    text = element_text(size = MY_FONT_SIZE, family = MY_FONT_FAMILY),
    legend.key = element_rect(color = NA),
    legend.key.width = unit(1, "lines"),
    legend.margin = margin(0.15, 0.5, 0.15, 0.3, "lines"),
    #legend.background = element_rect(linewidth = 0.5, color = "gray65"),
    legend.text = element_text(size = MY_FONT_SIZE_LEGEND))
}

my_labeller_np <- function(df, x, y) {
  if (missing(x)) x <- (1 + length(levels(df$K)))/2
  
  df_label <- df %>%
    select(method, n, p, setting_id) %>%
    distinct() %>%
    group_by(method, n, p, setting_id) %>%
    summarise(label   = paste0("n = ", n, ", dim(X) = ", p),
              n_label = paste0("n = ", n),
              p_label = paste0("dim(X) = ", p), .groups = "keep") %>%
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

my_labeller_trees <- function(df, x, y) {
  if (missing(x)) x <- (1 + length(levels(df$K)))/2
  
  df_label <- df %>%
    select(method, nt, setting_id) %>%
    distinct() %>%
    group_by(method, nt, setting_id) %>%
    summarise(label   = paste0("Trees = ", nt),
              nt_label = paste0("Trees = ", nt),
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

df_plt <- df_ratio %>% filter(model_type == MODEL_TYPE)

### MODEL-SPECIFIC SETTINGS
my_ylim <- switch(MODEL_TYPE,
                  "vcm" = c(0.6, 3.15),
                  "hte" = c(0.85, 1.65))
my_breaks <- pretty(seq(min(my_ylim), max(my_ylim), length.out = 5))

#my_label_ypos <- switch(MODEL_TYPE, "vcm" = 1, "hte" = 1)
my_label_np_ypos <- switch(MODEL_TYPE, "vcm" = 0.75, "hte" = 0.90)
my_label_trees_ypos <- switch(MODEL_TYPE, "vcm" = 3.05, "hte" = 1.60)
my_model_colors <- switch(MODEL_TYPE, "vcm" = MY_COLORS[1:4], "hte" = MY_COLORS)

title_str <- "Forest fit time speedup factor: GRF-grad time/GRF-FPT time"
subtitle_str <- 
  switch(MODEL_TYPE,
         "vcm" = "Varying coefficient model (VCM)",
         "hte" = "Heterogeneous treatment effects (HTE)")
legend_str <- sprintf("%s\nSetting", toupper(MODEL_TYPE))

### MAKE PLOTS (line/point)
plt_line <- df_plt %>%
  ggplot(aes(x = K, y = ratio, group = setting_id, color = setting_id, shape = setting_id)) +
  geom_hline(aes(yintercept = hline), col = 'gray75', linewidth = 0.75, linetype = "solid") +
  labs(x = "Number of regressors (K)", y = "Speedup factor", color = legend_str, shape = legend_str) +
  ggtitle(title_str, subtitle = subtitle_str) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(nt ~ n + p) +
  scale_y_continuous(limits = my_ylim, breaks = my_breaks) +
  scale_color_manual(values = my_model_colors) + 
  my_theme() +
  my_labeller_np(df_plt, y = my_label_np_ypos) +
  my_labeller_trees(df_plt, y = my_label_trees_ypos)

### MAKE PLOTS (barplot)
df_plt_adjusted <- df_plt %>%
  mutate(ratio_from_one = ratio - 1) # calculate height from baseline of 1

plt_bar <- df_plt_adjusted %>%
  ggplot(aes(x = K, y = ratio_from_one, fill = setting_id)) + 
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_hline(yintercept = 0, col = 'gray75', linewidth = 0.75, linetype = "solid") +
  labs(x = "Regressor dimension (K)", y = "Speedup factor", fill = legend_str) +
  ggtitle(title_str, subtitle = subtitle_str) + 
  facet_grid(nt ~ n + p) +
  scale_y_continuous(
    limits = my_ylim - 1,
    labels = function(x) sprintf("%.1f", x + 1),
    breaks = my_breaks - 1) +
  scale_fill_manual(values = my_model_colors) + 
  my_theme() +
  my_labeller_np(df_plt, y = my_label_np_ypos - 1) +
  my_labeller_trees(df_plt, y = my_label_trees_ypos - 1)

#filename_plt_line <- sprintf("figures/forest-timing-ratio-%s-line.png", MODEL_TYPE)
filename_plt_bar <- sprintf("figures/forest-timing-ratio-%s-bar.png", MODEL_TYPE)

#ggsave(filename_plt_line, plot = plt_line, width = 6, height = 7)
ggsave(filename_plt_bar, plot = plt_bar, width = 7, height = 4.5)


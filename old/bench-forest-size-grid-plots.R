##############################################################################################################
# forest-grid-VCM-plots.R
#
#
#
# Generated on 2025-01-14 14:47 EST
##############################################################################################################
library(tidyverse)
library(ggh4x)
#load("data/bench-forest-size-grid-vcm-data-20250115-2057.rdata")
load("data/forest-grid-vcm-sim-20250115-2057.rdata")
# Glance at the unique values
df_main %>%
  select(method, K, p, n, nt, mns, rep) %>%
  map(unique)
df_main %>%
  select(method, K, p, n, nt, mns, rep) %>%
  map(table)

#grouping_vars <- c("method", "K", "p", "n", "nt", "mns", "rep")
grouping_vars <- c("method", "K", "p", "n", "nt", "mns")

df_summary_nsplits <- df_nsplits %>%
  group_by(across(all_of(grouping_vars))) %>% 
  summarise(nsplits = mean(nsplits), .groups = "keep")

df_summary <- df_main %>%
  group_by(across(all_of(grouping_vars))) %>% 
  summarise(time = median(time), 
            amse = mean(amse), .groups = "keep") %>%
  inner_join(df_summary_nsplits, 
             by = join_by(!!!syms(grouping_vars)))

df_ratio <- df_summary %>%
  group_by(across(all_of(setdiff(grouping_vars, "method")))) %>%
  summarise(
    across(c(time, amse, nsplits),
           list(FPT1 = ~.[method == "FPT1"] / .[method == "grad"],
                FPT2 = ~.[method == "FPT2"] / .[method == "grad"])), 
    .groups = "drop") %>%
  pivot_longer(cols = contains("_"),
               names_to = c("metric", "method"),
               names_pattern = "(.*)_(.*)",
               values_to = "value") %>%
  pivot_wider(names_from = metric,
              values_from = value) %>%
  mutate(hline = 1) # for plotting later

# # Compute relative benchmark times/relative AMSE
# df_ratio_main <- df_summary %>%
#   group_by(K, p, n, nt, mns, rep) %>% {
#     time_frame <- reframe(.,
#                           "FPT1" = time[method == "FPT1"]/time[method == "grad"],
#                           "FPT2" = time[method == "FPT2"]/time[method == "grad"]) %>%
#       gather(method, time, "FPT1":"FPT2")
#     amse_frame <- reframe(.,
#                           "FPT1" = amse[method == "FPT1"]/amse[method == "grad"],
#                           "FPT2" = amse[method == "FPT2"]/amse[method == "grad"]) %>%
#       gather(method, amse, "FPT1":"FPT2")
#     by <- join_by(K, p, n, nt, mns, rep, method)
#     inner_join(time_frame, amse_frame, by = by)
#   } %>% mutate(hline = 1) # for plotting later

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

my_theme <- function(legend_pos = NULL, legend_loc = "inside") {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text = element_text(family = MY_FONT_FAMILY),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.major.x = element_line(color = "gray95"),
    legend.position = legend_loc,
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
my_common_elements <- function(xlab = "Regressors (K)", ylab, 
                               plottype = c("point", "boxplot"), ...) {
  res <- list(xlab(xlab),
              ylab(ylab),
              ggh4x::facet_nested("Minimum node size (target)" + mns ~ "Number of trees" + nt,
                                  nest_line = element_line(colour = "gray50")))
              # facet_grid(vars(nt), vars(mns)))
  if (plottype == "point") {
    res <- c(list(geom_line(linewidth = 0.5),
                  geom_point(size = 2.5, stroke = 0.75)), res)
  } else if (plottype == "boxplot") {
    res <- c(list(geom_boxplot(alpha = 0.2), res))
  } else {
    stop("Something went wrong in my_common_elements(). Check 'plottype'.")
  }
  return (res)
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
# Optional if plots are too ugly
nt_keep <- unique(df_ratio$nt) 
mns_keep <- unique(df_ratio$mns)
df_ratio_plot <- df_ratio %>%
  filter(nt %in% nt_keep, mns %in% mns_keep) %>%
  filter(mns != 640)


plottype <- if ("rep" %in% colnames(df_ratio)) "boxplot" else "point"
my_group <- if (plottype == "boxplot") NULL else quo(method)

plt_nsplits <- df_ratio_plot %>% 
  ggplot(aes(x = K, y = 1/nsplits, shape = method, fill = method, color = method, group = !!my_group)) +
  my_ratio_elements(plottype = plottype, ylab = "Split ratio") + 
  ggtitle("Relative tree complexity: (avg splits GRF-grad) / (avg splits GRF-fpt)") +
  my_scale_y(expand = c(0.05, 0.05)) + 
  #my_labeller(df_ratio, y = 0.25) +
  my_theme(legend_loc = "bottom")
plt_nsplits

plt_amse <- df_ratio_plot %>% 
  ggplot(aes(x = K, y = 1/amse, shape = method, fill = method, color = method, group = !!my_group)) +
  my_ratio_elements(plottype = plottype, ylab = "avgMSE factor") + 
  ggtitle("Relative avgMSE factor (median GRF-grad avgMSE) / (median GRF-fpt avgMSE)") +
  #my_scale_y(expand = c(0.05, 0.05)) + 
  #my_labeller(df_ratio, y = 0.25) +
  my_theme(legend_loc = "bottom")
plt_amse

plt_time <- df_ratio_plot %>% 
  ggplot(aes(x = K, y = 1/time, shape = method, fill = method, color = method, group = !!my_group)) +
  my_ratio_elements(plottype = plottype, ylab = "Median speedup factor") + 
  ggtitle("Timing benchmarks (forests): Fit time speedup factor (median GRF-grad time) / (median GRF-fpt time)") +
  #my_scale_y(expand = c(0.05, 0.05)) + 
  #my_labeller(df_ratio, y = 0.25) +
  my_theme(legend_loc = "bottom")
plt_time

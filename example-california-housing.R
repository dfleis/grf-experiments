##############################################################################################################
# example-california-housing.R
# Real data example of GRF estimates using the California housing data.
#
# NOTES (Data):
#   * See https://www.dcc.fc.up.pt/~ltorgo/Regression/cal_housing.html
#   * Originally (?) appearing in the paper 
#       Pace, R.K., and Barry, R. (1997). Sparse Spatial Autoregressions. 
#         Statistics and Probability Letters. 33(3), 291-297.
#     For the paper, see
#       https://doi.org/10.1016/S0167-7152(96)00140-X
#       https://www.sciencedirect.com/science/article/pii/S016771529600140X
#   * High-level summary (from ageron/handson-ml/tree/master/datasets/housing)
#       The dataset was derived from 1990 California census data. Each row ("observation") 
#       corresponds to the census data associated with a census "block group". A census
#       block group is the smallest geographical unit for which the U.S. Census Bureau
#       publishes sample data (typically has a population of 600 to 3000 people).
#
# NOTES (Model):
#   * Varying coefficient model, e.g. Hastie and Tibshirani (1993)
#     http://www.jstor.org/stable/2345993 of the form
#       E[Y_i | X_i = x] = \nu^*(x) + W_{i,1} \theta^*_1(x) + ... + W_{i,K} \theta^*_K(x), 
#                        = \nu^*(x) + W_i^\top \theta^*(x),
#     for primary regressors W_i = (W_{i,1}, ..., W_{i,K})^\top \in \mathbb R^K, 
#     auxiliary/mediating covariates X_i \in \mathbb R^p, outcome/response Y_i \in \mathbb R. 
#     The function \nu^*(x) is a nuisance intercept function and our target is estimates of
#     the local effects \theta^*(x) = (\theta^*_1(x), ..., \theta^*_K(x))^\top).
#   * The model we use for this example is a spatially-varying coefficient model for median
#     house value:
#       X = (latitude, longitude),
#       W = (... regressors/features/attributes ...),
#       Y = (log) median house value
#
# Generated 2025-Jan-06 23:05 EST
###################################################################################################################
set.seed(1)
library(tidyverse)
library(GGally)
library(bench)
library(grf)
library(sf)
library(tigris) # state boundaries
options(tigris_use_cache = TRUE)

PLOT_PATH <- "figures/california-housing"

do_sprintf <- function(.fmt, .args) {
  do.call(sprintf, c(.fmt, as.list(.args)))
}
my_ggsave <- function(plot, filename, ext = c("pdf", "png"), path = NULL,
                      width, height, units = "in", bg = "white", 
                      create.dir = TRUE, ...) {
  if (is.null(path)) path <- getwd()
  fnms <- paste(filename, ext, sep = ".")
  lapply(fnms, function(ff) {
    cat(format(Sys.time()), sprintf("ggsave %s/%s\n", path, ff))
    ggsave(filename = ff, plot = plot, path = path,
           width = width, height = height, units = units, 
           bg = bg, create.dir = create.dir, ...)
  })
  invisible()
}

#-------------------------------------------------------
#----- PREPARE DATA
#----- Source: https://www.dcc.fc.up.pt/~ltorgo/Regression/cal_housing.tgz
#-------------------------------------------------------
# Load data
data.colnames <- file("data/CaliforniaHousing/cal_housing.domain") %>%
  readLines() %>%
  sub(pattern = ":.*", replacement = "", x = .)
data <- read.delim("data/CaliforniaHousing/cal_housing.data", 
                   sep = ",", header = F, col.names = data.colnames)

# Rename columns
data <- data %>%
  select(lon = longitude,
         lat = latitude,
         med_val = medianHouseValue,
         med_age = housingMedianAge,
         tot_rms = totalRooms,
         tot_bds = totalBedrooms,
         pop     = population,
         hholds  = households,
         med_inc = medianIncome)

# Separate the columns/data attributes according to their role in the model, i.e.
# covariates X_i, outcome Y_i, regressors W_i
X <- data %>% select(lat, lon) %>% as.matrix()
Y <- data %>% transmute(med_val = log(med_val)) %>% as.matrix()
W <- data %>%
  transmute(med_age = med_age,
            tot_rms = log(tot_rms),
            tot_bds = log(tot_bds),
            pop     = log(pop),
            hholds  = log(hholds),
            med_inc = log(med_inc)) %>%
  as.matrix()

W_labs <- c(
  med_age = "Median Housing Age",
  tot_rms = "log Census Block Rooms",
  tot_bds = "log Census Block Bedrooms",
  pop     = "log Population",
  hholds  = "log Households",
  med_inc = "log Median Income")

# #-------------------------------------------------------
# #----- STUMP COMPARISONS
# #-------------------------------------------------------
# seed <- 1
# methods <- c("grad", "fpt1", "fpt2")
# num.trees <- 2000 # number of stumps
# min.node.size <- nrow(X)/2 - 1 # enforce stumps
# 
# forests <- sapply(methods, function(method) {
#   forest <- lm_forest(X = X, Y = Y, W = W, 
#                       num.trees = num.trees, 
#                       min.node.size = min.node.size, 
#                       honesty = F,
#                       method = method, 
#                       seed = seed)
# }, USE.NAMES = TRUE, simplify = F)
# 
# df_stumps <- lapply(forests, function(ff) {
#   splitvar <- sapply(ff$`_split_vars`, `[[`, 1) + 1 # indexing starts at 0
#   splitval <- sapply(ff$`_split_values`, `[[`, 1)
#   data.frame(method = ff$method, splitvar = splitvar, splitval = splitval)
# }) %>% 
#   bind_rows() %>%
#   mutate(method = factor(method, levels = methods)) %>%
#   mutate(method = recode(method, fpt1 = "FPT1", fpt2 = "FPT2"))
# 
# splitvar_labels <- c("1" = "Latitude", "2" = "Longitude") 
# 
# plt_stumps <- df_stumps %>%
#   ggplot(aes(x = method, y = splitval, fill = method, color = method)) + 
#   geom_boxplot(alpha = 0.25, linewidth = 0.75, width = 0.5, 
#                outlier.alpha = 0.1, outlier.shape = 16) + 
#   facet_grid(rows = vars(splitvar), scales = "free", 
#              labeller = labeller(splitvar = splitvar_labels)) + 
#   labs(x = NULL, y = "Split Value",
#        title = "California housing data: Distribution of root node splits",
#        subtitle = sprintf("%s root node splits", num.trees),
#        color = "Method", fill = "Method") +
#   theme_minimal() + 
#   theme(
#     panel.border = element_rect(color = "black", fill = NA), 
#     axis.ticks.length = unit(-0.2, "cm"),
#     axis.ticks = element_line(),
#     axis.text = element_text(family = "mono", size = 12),
#     axis.line = element_line(),
#     strip.text = element_text(size = 12),
#     legend.position = "inside",
#     legend.position.inside = c(0.925, 0.6),
#     legend.title.align = 0.5,
#     legend.background = element_rect(fill = "white", color = "black"),
#     legend.text = element_text(family = "mono"))
# 
# my_ggsave(plt_stumps, "grf-california-housing-split-stability", 
#           path = PLOT_PATH, width = 7, height = 4)

#-------------------------------------------------------
#----- FIT FORESTS
#-------------------------------------------------------
seed <- 1
num.threads <- NULL
methods <- c("grad", "fpt1", "fpt2")
num.trees <- 2000
min.node.size <- 5

forests <- sapply(methods, function(method) {
  gc()
  tm <- bench::system_time(
    forest <- lm_forest(X = X, Y = Y, W = W, 
                        num.trees = num.trees, 
                        min.node.size = min.node.size, 
                        method = method, 
                        seed = seed, 
                        num.threads = num.threads)
  )
  return (list(forest = forest, time = tm))
}, USE.NAMES = TRUE, simplify = F)

#-------------------------------------------------------
#----- PREDICT NEW VALUES
#----- This script is exceptionally ugly, please 
#----- be gentle
#-------------------------------------------------------
### Part 1: Create grid of new X covariates (lat/lon values)
# define bounding box of the target region
bbox <- list( 
  lon_min = min(X[,"lon"]), # western extent
  lon_max = max(X[,"lon"]), # eastern extent
  lat_min = min(X[,"lat"]), # southern extent
  lat_max = max(X[,"lat"])) # northern extent

n_lat <- n_lon <- 200
coord_grid <- expand.grid(lat = seq(bbox$lat_min, bbox$lat_max, length.out = n_lat), 
                          lon = seq(bbox$lon_min, bbox$lon_max, length.out = n_lon))

# get California polygon
ca <- states() %>% 
  filter(STUSPS == "CA") %>%
  st_transform(4326)  # Transform to WGS84

# convert grid points to sf object
points_sf <- st_as_sf(coord_grid, coords = c("lon", "lat"), crs = 4326)

# keep only points within California
ca_points <- points_sf[ca,]

# convert back to data frame
ca_coords <- data.frame(
  lon = st_coordinates(ca_points)[,1],
  lat = st_coordinates(ca_points)[,2]
)

Xnew_grid <- cbind(lat = ca_coords$lat, lon = ca_coords$lon)

########################################################
for (MAKE_INSAMPLE_PREDS in c(FALSE, TRUE)) {
# MAKE_INSAMPLE_PREDS <- FALSE
# MAKE_INSAMPLE_PREDS <- TRUE
cat(sprintf("%s :: Making predictions :: in-sample predictions = %s\n", 
            format(Sys.time()), MAKE_INSAMPLE_PREDS))
########################################################

### Part 2: Make predictions at X values
df_preds_list <- sapply(forests, function(ff) {
  if (MAKE_INSAMPLE_PREDS) {
    Xnew <- X  
    preds <- predict(ff$forest, estimate.variance = TRUE)
    pp <- preds$predictions[,,1] 
    #pp <- preds$variance.estimates
  } else {
    Xnew <- Xnew_grid
    preds <- predict(ff$forest, newdata = Xnew, estimate.variance = TRUE)
    pp <- preds$predictions[,,1]
    #pp <- preds$variance.estimates
  }
  return (data.frame(method = ff$forest$method,
                     lat = Xnew[,"lat"],      
                     lon = Xnew[,"lon"], 
                     pp))
}, USE.NAMES = TRUE, simplify = FALSE)


df_preds_long <- df_preds_list %>%
  bind_rows() %>%
  pivot_longer( # covert wide-formatted data to long-formatted data
    cols = any_of(colnames(W)),
    names_to = "effect",
    values_to = "value") %>%
  mutate(method = factor(method, levels = methods))

#-------------------------------------------------------
#----- VISUALIZE RESULTS
#-------------------------------------------------------
# get state boundaries
states_map <- map_data("state")
ca_map <- subset(states_map, region == "california")
major_cities <- data.frame(
  name = c("San Francisco", "Los Angeles", "Sacramento", "San Diego"),
  lon = c(-122.4194, -118.2426, -121.4944, -117.1611),
  lat = c(37.7749, 34.0549, 38.5781, 32.7157),
  nudge_x = c(2, -1.3, 1.3, 1.2),
  nudge_y = c(-0.2, -0.4, 0.4, 0.4)
)

########################################################
for (METHOD_TO_PLOT in methods) {
# METHOD_TO_PLOT <- "grad"
# METHOD_TO_PLOT <- "fpt1"
# METHOD_TO_PLOT <- "fpt2"
cat(sprintf("%s :: Plotting method = %s\n", 
            format(Sys.time()), METHOD_TO_PLOT))
########################################################

title <- sprintf("California housing data: GRF-%s", 
                 ifelse(METHOD_TO_PLOT == "grad", METHOD_TO_PLOT, toupper(METHOD_TO_PLOT)))
subtitle <- "Spatially-varying effects on log median house values"
if (MAKE_INSAMPLE_PREDS) subtitle <- paste0(subtitle, " (in-sample predictions)")

plot_xlim <- range(df_preds_long$lon)
plot_ylim <- range(df_preds_long$lat)
plot_colorlim <- range(df_preds_long$value)
plot_colorlim <- c(-max(abs(plot_colorlim)), max(abs(plot_colorlim)))
df_plot <- df_preds_long %>%
  filter(method == METHOD_TO_PLOT)# %>%
  #filter(effect %in% c("tot_rooms", "tot_beds", "med_age", "med_income"))

point_shape <- ifelse(MAKE_INSAMPLE_PREDS, 16, 15)
point_size <- ifelse(MAKE_INSAMPLE_PREDS, 0.4, 0.2)
 
# POINT_COLOR_LOW <- "#00BFC4"
# POINT_COLOR_HIGH <- "#F8766D"

color_transition_factor <- 1.5
color_pal <- function(n = 100) {
  # apply a tanh transformation to make the color gradient change more rapidly near zero
  vals <- seq(-1, 1, length.out = n)
  vals_transformed <- tanh(color_transition_factor * vals)  # Apply tanh transformation
  # normalize to [0,1]
  vals_normalized <- (vals_transformed - min(vals_transformed)) / 
    (max(vals_transformed) - min(vals_transformed))
  # create color palette
  #colorRamp(colors = c(POINT_COLOR_LOW, "white", POINT_COLOR_HIGH))(vals_normalized)
  colorRamp(colors = rev(pals::brewer.spectral(n)))(vals_normalized)
  colorRamp(colors = rev(pals::brewer.rdylbu(n)))(vals_normalized)
}

plt_preds <- ggplot() +
  geom_point(data = df_plot, size = point_size, shape = point_shape,
             aes(x = lon, y = lat, color = value, fill = value)) + 
  facet_wrap(vars(effect), labeller = labeller(effect = W_labs)) + 
  labs(x = "", y = "", color = "Effect on\nlog value", fill = "Effect on\nlog value",
       title = title, subtitle = subtitle) + 
  geom_point(data = major_cities, aes(x = lon, y = lat), color = "black", size = 1) +
  geom_sf(data = ca, fill = NA) +
  coord_sf(xlim = plot_xlim, ylim = plot_ylim) + 
  geom_label(data = major_cities, aes(x = lon, y = lat, label = name),
             label.size = 0, label.r = unit(0, "lines"),
             fill = "white", color = "black", size = 2.5, 
             nudge_x = major_cities$nudge_x,
             nudge_y = major_cities$nudge_y) + 
  # scale_color_gradient2(low = POINT_COLOR_LOW, high = POINT_COLOR_HIGH, mid = "white", 
  #                       midpoint = 0, limits = plot_colorlim, trans = color_scale_transform) +
  # scale_fill_gradient2(low = POINT_COLOR_LOW, high = POINT_COLOR_HIGH, mid = "white", 
  #                      midpoint = 0, limits = plot_colorlim, trans = color_scale_transform) +
  scale_color_gradientn(colours = rgb(color_pal(), maxColorValue = 255),
                        limits = plot_colorlim) +
  scale_fill_gradientn(colours = rgb(color_pal(), maxColorValue = 255), 
                       limits = plot_colorlim) +
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", margin = margin(b = 1)),
        axis.text = element_blank(), 
        legend.position = "inside",
        #legend.title.align = 0.5,
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position.inside = c(0.03, 0.63),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(1, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        panel.spacing.y = unit(0, "lines"),
        plot.margin = unit(c(0, 0, -1.5, 0), "lines"))  # top, right, bottom, left

filename <- do_sprintf("grf-california-housing-effects%s-%s", 
                       c(ifelse(MAKE_INSAMPLE_PREDS, "-insample", ""), METHOD_TO_PLOT))
my_ggsave(plt_preds, filename, path = PLOT_PATH, width = 7, height = 6.5)
} # end METHOD_TO_PLOT loop
} # end MAKE_INSAMPLE_PREDS loop

#-------------------------------------------------------
#----- EXTRA FIGURES: DISTRIBUTION OF THE TRAINING DATA 
#-------------------------------------------------------
set.seed(1)
W.plot.idx <- sample(nrow(W), 2000, replace = F)
plt_W <- data.frame(W[W.plot.idx,]) %>%
  rename(!!!setNames(names(W_labs), W_labs)) %>%
  ggpairs(title = "California housing data: Regressor distribution",
          labeller = label_wrap_gen(14),
          lower = list(continuous = wrap("points", alpha = 0.15, size = 0.1)),
          diag = list(continuous = wrap("barDiag", bins = 12, 
                                        fill = "gray95", color = "black")), 
          upper = list(continuous = wrap("cor", align_percent = 0.5,
                                         size = 3, stars = F))) +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA), 
        strip.text = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.ticks.length = unit(-0.1, "cm"),
        axis.ticks = element_line())

my_ggsave(plt_W, "grf-california-housing-regressor-distribution", 
          path = PLOT_PATH, width = 6, height = 6.2)

###################################################################
sapply(forests, `[[`, "time")


# #-------------------------------------------------------
# #----- TIME TESTS
# #-------------------------------------------------------
# seed <- 1
# num.threads <- NULL
# methods <- c("grad", "fpt1", "fpt2")
# num.trees <- 2000
# min.node.size <- 5
# 
# library(bench)
# 
# t0 <- Sys.time()
# bp <- bench::press(
#   method = methods,
#   {
#     bench::mark(
#       lm_forest(X = X, Y = Y, W = W, 
#                 num.trees = num.trees, 
#                 min.node.size = min.node.size, 
#                 method = method, 
#                 seed = seed, 
#                 num.threads = num.threads),
#       iterations = 10,
#       filter_gc = FALSE # large forests appear to almost always trigger GC
#     )
#   }
# )
# t1 <- Sys.time()
# difftime(t1, t0)
# 
# bp
# knitr::kable(bp %>% select(method, median))


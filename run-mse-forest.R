suppressMessages(library(tidyr))
suppressMessages(library(purrr))
suppressMessages(library(dplyr))
library(grf)
source("utils/utils-data.R")
source("utils/utils-forest.R")

#--------------------------------------------------
#----- Settings, parameters
#--------------------------------------------------
seed <- 1
NUM_THREADS <- NULL 
FILENAME_HEAD <- "forest-mse"

mse_nreps <- 50

#----- Data-generating parameters
args <- commandArgs(TRUE)
model_type <- args[1] # "vcm" or "hte"
setting_id <- args[2] # "1", "2", "3", "4", "5"

stopifnot(model_type %in% c("vcm", "hte"))
stopifnot(setting_id %in% c("1", "2", "3", "4", "5"))
settings <- validate_setting(model_type = model_type, setting_id = setting_id)

#----- Data params
Kvals <- c(4, 256)
pvals <- 5
nvals <- as.integer(c(1e4, 2e4, 1e5))
n_test <- 5000

#----- GRF parameters
num_trees <- c(100)
methods <- c("grad", "fpt1", "fpt2")

# Global GRF arguments
grf_args_global <- list(
  min.node.size = 5,
  sample.fraction = 0.5,
  honesty = TRUE,
  honesty.fraction = 0.5,
  honesty.prune.leaves = TRUE,
  alpha = 0.05,
  ci.group.size = 1,
  compute.oob.predictions = FALSE,
  num.threads = NUM_THREADS,
  seed = seed
)

# Note on ci.group.size = 1:
#
# Each subsample is used to train ci.group.size trees. The conventional
# presentation of the GRF algorithm corresponds to ci.group.size = 1
# because every tree uses its own subsample. However, for CI estimation,
# GRF will compute trees in small batches were an individual subsample 
# is used to train ci.group.size trees.
#
# This mechanism affects the random number generator used in the
# C++ backend to do the subsampling in a way that makes it difficult
# to guarantee that two forests will use the same collection of 
# subsamples for training their trees. Specifically, in order for the
# drawn subsamples to be identical across the different methods 
# ("grad", "fpt1", "fpt2"), we disable CI batch sampling by setting
# ci.group.size = 1. This is fine for this simulation because we're 
# studying the forest fit times and not CI coverage.

#--------------------------------------------------
#----- Simulations
#--------------------------------------------------
set.seed(seed)
par_grid <- expand.grid(K = Kvals, p = pvals, n = nvals, nt = num_trees)
FUN_grf <- get_FUN_grf(model_type) 

t0 <- Sys.time()
sim_list <- vector(mode = "list", length = nrow(par_grid))
for (i in 1:nrow(par_grid)) {
  pars <- par_grid[i,]
  K <- pars$K
  p <- pars$p
  n <- pars$n
  nt <- pars$nt
  
  cat(format(Sys.time()), "::", i, "of", nrow(par_grid), "::",
      paste(model_type, setting_id), "::", FUN_grf$name, "::",
      paste(names(pars), pars, sep = " = ", collapse = ", "), "\n")
  
  df_list <- lapply(1:mse_nreps, function(repid) {
    cat(format(Sys.time()), ":: repid =", repid, "\n")
    # Generate data
    data_train <- generate_data(
      n = n, p = p, K = K, 
      model_type = settings$model_type,
      effect_type = settings$effect_type, 
      prob_type = settings$prob_type
    )
    X_test <- generate_X(n_test, p)
    thetaX_test <- data_train$theta(X_test)

    # Common centering step as per the recommendations of grf::lm_forest and grf::multi_arm_causal_forest
    centered <- make_centering_data(
      data = data_train, 
      model_type = model_type, 
      center_data = TRUE, 
      grf_args = modifyList(grf_args_global, list(num.trees = nt))
    )
    
    args_grf <- modifyList(grf_args_global, c(centered, list(num.trees = nt)))
    
    make_args <- function(mm) modifyList(args_grf, list(method = mm))

    # Fit forest and make predictions
    pred_stats <- sapply(methods, function(method) {
      # Fit forest (Stage I)
      forest <- do.call(FUN_grf$FUN, make_args(method))
      
      # Make predictions (Stage II)
      preds <- predict(forest, newdata = X_test, drop = TRUE)
      est <- preds$predictions
      
      colnames(est) <- colnames(thetaX_test)
      err <- est - thetaX_test
      dim_mse <- colMeans(err^2) # component-wise MSE
      avg_mse <- mean(dim_mse)   # average MSE
      nsplits <- (sapply(forest$`_split_vars`, length) - 1)/2 
      
      df_dim_mse <- data.frame(effect = names(dim_mse), dim_mse = dim_mse, row.names = NULL)
      df_nsplits <- data.frame(tree = 1:length(nsplits), nsplits = nsplits)
      
      list(avg_mse = avg_mse, dim_mse = df_dim_mse, nsplits = df_nsplits)
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    pred_stats_list <- purrr::list_transpose(pred_stats)  
    
    data.frame(method = methods) %>%
      mutate(pars, 
             rep = repid,
             model_type = settings$model_type,
             setting_id = settings$setting_id,
             effect_type = settings$effect_type,
             prob_type = settings$prob_type,
             FUN = FUN_grf$name,
             avg_mse = unname(pred_stats_list$avg_mse),
             dim_mse = pred_stats_list$dim_mse,
             nsplits = pred_stats_list$nsplits)
  })
  
  sim_list[[i]] <- bind_rows(df_list)
}
t1 <- Sys.time()

#--------------------------------------------------
#----- Post-processing
#--------------------------------------------------
df_sim <- bind_rows(sim_list) %>%
  mutate(method = factor(method, levels = methods))# %>%
  #mutate(method = recode(method, fpt2 = "FPT"))

df_avg_mse <- df_sim %>% select(-dim_mse, -nsplits)
df_dim_mse <- df_sim %>% 
  select(-avg_mse, -nsplits) %>%
  unnest_longer(dim_mse) %>%
  unnest(dim_mse)
df_nsplits <- df_sim %>%
  select(-avg_mse, -dim_mse) %>%
  unnest_longer(nsplits) %>%
  unnest(nsplits)

#--------------------------------------------------
#----- Save data
#--------------------------------------------------

datetime <- format(Sys.time(), format = "%Y%m%d-%H%M")

make_filename <- function(label) {
  file_dir <- sprintf("data/%s", gsub("/$", "", FILENAME_HEAD))
  dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)
  sprintf("%s/%s-%s-%s-%s-%s.csv", file_dir, FILENAME_HEAD, label, model_type, setting_id, datetime)
}

write.csv(df_avg_mse, file = make_filename("avg_mse"), row.names = FALSE)
write.csv(df_dim_mse, file = make_filename("dim_mse"), row.names = FALSE)
write.csv(df_nsplits, file = make_filename("nsplits"), row.names = FALSE)

cat("Total elapsed time:", format(difftime(t1, t0), digits = 6), "\n")

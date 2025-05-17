suppressMessages(library(tidyverse))
library(bench)

source("utils/utils-data.R")
source("utils/utils-forest.R")

#-------------------------------------------------------
#----- Post-processing tools
#-------------------------------------------------------
process_bench_forest <- function(bench_data, data_pars) {
  drop_bench_cols <- c("expression", "min", "median", "itr/sec", 
                       "mem_alloc", 
                       "gc/sec", "n_itr", "n_gc", 
                       "total_time", "result", "memory")
  
  bench_data_small <- bench_data %>%
    select(-all_of(drop_bench_cols))
  
  splits_list <- lapply(bench_data$result, function(forest) {
    splits <- (sapply(forest$`_split_vars`, length) - 1)/2 # number of splits per tree
    data.frame(tree = 1:length(splits), splits = splits)
  })
  
  df_splits <- bench_data_small %>%
    mutate(splits = splits_list) %>%
    select(-time, -gc) %>%
    unnest_longer(splits) %>%
    unnest_wider(splits) %>%
    mutate(data_pars)
  
  df_time <- bench_data_small %>%
    unnest_longer(c(gc, time), indices_include = T) %>%
    unnest_wider(gc) %>%
    rename(gc0 = level0, gc1 = level1, gc2 = level2, iter = time_id) %>%
    mutate(gc   = gc0 + gc1 + gc2 > 0,
           time = as.numeric(time)) %>%
    select(method, rep, iter, gc, gc0, gc1, gc2, time) %>%
    mutate(data_pars)
  
  list(splits = df_splits, time = df_time)
}

#-------------------------------------------------------
#----- Benchmarking tools
#-------------------------------------------------------
bench_forest_pars <- function(methods, model_type, 
                              data_pars, center_data = FALSE, 
                              grf_args = NULL,
                              nrep, niter, 
                              seed = NULL, .quiet = TRUE) {
  # model_type: "vcm", "hte"
  # methods: "grad", "fpt1", "fpt2"
  # data_pars: setting_id, K, p, n
  # nrep: bench replications (see ?bench::mark and ?bench::press)
  # niter: bench iterations (see ?bench::mark and ?bench::press)
  grf_FUN <- get_FUN_grf(model_type) # vcm: lm_forest, hte: multi_arm_causal_forest
  
  res <- bench::press(
    rep = 1:nrep, # dummy variable to do replications over new data samples
    method = methods,
    {
      data_args <- data_pars
      data_args$model_type <- model_type
      data_args$seed <- ifelse(is.null(seed), rep - 1, seed + rep - 1)
      
      data <- do.call(generate_data, data_args)
      centering_data <- make_centering_data(data = data, 
                                            model_type = model_type, 
                                            center_data = center_data, 
                                            grf_args = grf_args)
      
      grf_bench_args <- modifyList(grf_args, centering_data)
      grf_bench_args$method <- method
      
      bench::mark(
        do.call(grf_FUN$FUN, grf_bench_args), 
        iterations = niter,
        min_time = Inf,
        filter_gc = FALSE,
        memory = FALSE
      )
    },
    .quiet = .quiet
  )
  
  pars <- bind_cols(model_type = model_type, 
                    num.trees = grf_args$num.trees, 
                    data_pars)
  return (process_bench_forest(res, pars))
}

#--------------------------------------------------
#----- Core bench function
#--------------------------------------------------
bench_forest <- function(methods, model_type, 
                         setting_id, numtreevals, Kvals, pvals, nvals, 
                         center_data = FALSE, grf_args = NULL, 
                         nrep, niter, seed = NULL,
                         path = "data", filename_head = "bench-forest", 
                         .add_col = NULL,
                         .quiet = TRUE) {
  setting_id <- validate_setting(model_type, setting_id)$setting_id
  
  pars_grid <- expand.grid(
    K = sort(Kvals, decreasing = T),
    p = sort(pvals, decreasing = T), 
    n = nvals,#sort(nvals, decreasing = T),
    num.trees = sort(numtreevals, decreasing = T),
    setting_id = setting_id
  )
  
  file_dir <- sprintf("%s/%s-%s-%s", path, filename_head, nrep, niter)
  dir.create(file_dir, recursive = TRUE, showWarnings = FALSE)

  t0 <- Sys.time()
  for (i in 1:nrow(pars_grid)) {
    
    pars <- pars_grid[i,]
    pars_str <- paste(c(model_type, setting_id, as.integer(pars$num.trees), as.integer(pars$K), 
                        as.integer(pars$p), as.integer(pars$n)), collapse = "-")
    # file_dir/filename_head-[NREPS]-[NITERS]-{VARLABEL}-[MODEL]-[SETTING]-[NUM.TREES]-[K]-[p]-[n].csv
    make_filename <- function(LABEL) {
      sprintf("%s/%s-%s-%s-%s-%s.csv", file_dir, 
              filename_head, nrep, niter, LABEL, pars_str)
    }
    
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "::", 
        i, "of", nrow(pars_grid), "::", make_filename("{VARLABEL}"), "\n")
    gc(); Sys.sleep(1)
    
    res <- bench_forest_pars(
      methods = methods, 
      model_type = model_type, 
      data_pars = modifyList(pars, list(num.trees = NULL)),
      center_data = center_data,
      grf_args = modifyList(grf_args, list(num.trees = pars$num.trees)),
      nrep = nrep, 
      niter = niter, 
      seed = seed,
      .quiet = .quiet
    )
    
    if (!is.null(.add_col)) {
      res$splits <- res$splits %>% bind_cols(.add_col)
      res$time <- res$time %>% bind_cols(.add_col)
    }
    
    write.csv(res$splits, file = make_filename("SPLITS"), row.names = FALSE)
    write.csv(res$time, file = make_filename("TIME"), row.names = FALSE)
  }
  
  cat(sprintf("Benchmark complete for model_type = %s with\n", model_type))
  for (j in 1:ncol(pars_grid)) {
    cat("  ", colnames(pars_grid)[j], ":", unique(as.character(pars_grid[,j])), "\n")
  }
  cat(sprintf("Benchmark replications nrep = %s\n", nrep))
  cat(sprintf("Benchmark iterations niter = %s\n", niter))
  cat("Total elapsed time:", format(difftime(Sys.time(), t0)), "\n")
}

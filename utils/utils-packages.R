local({
  packages <- c(
    "devtools",  # installing `dfleis/rfg` and the FPT version of `dfleis/grf`
    "mvtnorm",   # multivariate normal
    "bench",     # benchmarking
    "tidyr",     # data manipulation
    "dplyr",     # data manipulation
    "purrr",     # data manipulation
    # Below are optional tools used for creating the figures
    "ggplot2",   # (optional) visualization
    "ggh4x",     # (optional) visualization
    "scales",    # (optional) visualization
    "pals"       # (optional) visualization
  )     
  
  install_and_load_packages <- function(pkgs) {
    for (pkg in pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing package:", pkg))
        install.packages(pkg, dependencies = TRUE)
      }
    }
    sapply(pkgs, require, character.only = TRUE)
  }
  
  install_and_load_packages(packages)  
})
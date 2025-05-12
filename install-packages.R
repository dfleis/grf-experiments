####################################################################################################
# install-packages.R
# Installs/builds and loads packages.
#
# WARNING: 
#   * Do not install/load `grf` through the standard `install.packages()` mechanism. Doing
#     so will install the CRAN version of `grf` which does not include the FPT method.
#
# Generated on 2025-05-11 20:24 EDT
####################################################################################################
#----- Load CRAN packages
cat("Loading utility packages:\n")
suppressMessages(source("utils/utils-packages.R", print.eval = TRUE))

#----- Load custom packages
if (!requireNamespace("rfg", quietly = TRUE)) {
  message(paste("Installing package `rfg` from `github.com/dfleis/rfg`"))
  devtools::install_github("dfleis/rfg")
}
success.rfg <- require(rfg, quietly = TRUE)
cat("Loaded `rfg`?", success.rfg, "\n")


if (!requireNamespace("grf", quietly = TRUE)) {
  message(paste("Installing package `grf` from `github.com/dfleis/grf`"))
  devtools::install_github("dfleis/grf", subdir = "r-package/grf")
}
success.grf <- require(grf, quietly = TRUE)
cat("Loaded `grf`?", success.grf, "\n")

if (!requireNamespace("grf", quietly = TRUE)) {
  message("Something went wrong installing the FPT `grf` library...")
} else {
  if (!isTRUE("method" %in% formalArgs(grf::lm_forest)) ||
      !isTRUE("method" %in% formalArgs(grf::multi_arm_causal_forest))) {
    stop("The installed version of `grf` does not contain the FPT implementation.",
         "\nPlease install the package from \"github.com/grf\" via",
         "\n\tdevtools::install_github(\"dfleis/grf\", subdir = \"r-package/grf\")")
  }
}
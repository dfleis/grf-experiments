# Build Instructions

## Unix-based systems

Installation of our fork of `grf` that implements the fixed-point tree algorithm can be done through the `remotes` or `devtools` libraries:
```R
# install.packages("remotes")
remotes::install_github("dfleis/grf", subdir = "r-package/grf")
```

# Windows

## Step 1: Install RTools

For R 4.1.x, use RTools41, for R 4.2.x, use RTools42, and so on.

## Step 2: Verify RTools installation

After installing RTools, ensure it is correctly configured. In R, run:
```R
Sys.getenv("PATH")
```

The output should include the path to RTools (e.g., for RTools44, something like ```C:\rtools44\usr\bin```). We can check that the compiler works by running:
```R
system("g++ --version")
```

## Step 3: Install packages to build grf
```R
install.packages("devtools") # install.packages("remotes")
install.packages("git2r")
```

## Step 4: Build grf

### Step 4a: Enable Windows developer mode

* In Windows, open Settings.
* Go to Update & Security > For Developers.
* Enable Developer Mode (this allows non-admin symlink creation).
* Run git config ```--global core.symlinks true```, or, in R, ```system("git config --global core.symlinks true")```.

### Step 4b: Clone the grf repo

```R
library(devtools)
library(git2r)
tmp_dir <- paste(getwd(), "tmp_dir", sep = "/")
dir.create(tmp_dir, recursive = TRUE)
#clone("https://github.com/grf-labs/grf.git", tmp_dir)
clone("https://github.com/dfleis/grf.git", tmp_dir)
```

### Step 4c: Build the package

```R
devtools::install(file.path(tmp_dir, "r-package", "grf"))
```

# Verify installation

```R
n <- 100
p <- 2
K <- 3

X <- matrix(rnorm(n * p), nrow = n)
W <- matrix(rnorm(n * K), nrow = n)
Y <- rnorm(n)
forest <- grf::lm_forest(X, Y, W, method = "fpt2")
forest$method
```

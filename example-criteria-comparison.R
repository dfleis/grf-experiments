####################################################################################################
# example-criteria-comparison.R
# Toy model to form a visual comparison of the target \Delta splitting criterion, the
# approximate \widetilde\Delta criterion specified by GRF-grad, and the approximate
# \tilde\Delta criteria specified by GRF-FPT under "step size" \eta = 1 and \eta = \sqrt{2}.
#
# Data-generating model (VCM):
#   Y_i = W_i^\top \theta^*(X_i) + \espilon_i
# where
#   * covariates X_i ~ Unif[0,1] (1-dimensional uniform),
#   * regressors W_i ~ N_2(0, I) (2-dimension Gaussian),
#   * noise \epsilon_i ~ N(0, 1),
#   * VCM effect functions \theta^*(x) = (\theta_1^*(x), \theta_2^*(x)) where
#       \theta_1^*(x) := \sin(2 \pi x),
#       \theta_2^*(x) := x.
#
####################################################################################################
library(ggplot2)
library(tidyr)
library(purrr)
library(dplyr)

#--------------------------------------------------
#----- Helper functions
#--------------------------------------------------
#--- VCM effect functions for K = 2
theta_FUN <- function(x) c(sin(2 * pi * x), x)

#--- VCM score function
psi_FUN <- function(theta, nu, yvec, wmat) {
  D <- cbind(1, wmat)
  params <- c(nu, theta)
  sweep(D, 1, yvec - D %*% params, FUN = "*")
}

#--- OLS coefficient calculator (used for \hat\theta_C)
lmcoef_FUN <- function(yvec, wmat) unname(coef(lm(yvec ~ wmat)))

#--- Child \hat\theta_C calculator
thetahatC_FUN <- function(C, yvec, wmat) {
  yy <- yvec[C]
  ww <- wmat[C,,drop=F]
  lmcoef_FUN(yy, ww)[-1]
}
#--- Child \tilde\theta_C calculator 
thetatildeC_FUN <- function(C, theta0, rho) theta0 + colMeans(rho[C,,drop=F])

#--- Gradient matrix estimator -A_P
A_FUN <- function(wmat) -crossprod(cbind(1, wmat))/nrow(wmat)

#--------------------------------------------------
#----- Generate data
#--------------------------------------------------
set.seed(1)
n <- 5000
K <- 2
sig.eps <- 0.25
sig.W <- 1

eta.new <- sqrt(3/2)

#--- "Selector" matrix such that \xi^\top (nu, \theta_1, \theta_2) = (\theta_1, \theta_2)
xi <- diag(c(0, rep(1, K)), nrow = K + 1, ncol = K + 1)[,-1]

#--- Covariates, regressors, noise, effects, outcomes
X <- matrix(sort(runif(n * 1)), nrow = n)
W <- matrix(rnorm(n * K, 0, sig.W), nrow = n)
eps <- rnorm(n, 0, sig.eps)
thetaX <- t(apply(X, 1, theta_FUN))
Y <- rowSums(thetaX * W) + eps

#--------------------------------------------------
#----- Compute criteria by definition
#--------------------------------------------------
lmcoefs.P <- lmcoef_FUN(Y, W)
thetahat.P <- lmcoefs.P[-1]
nuhat.P <- lmcoefs.P[1]

psi.P <- psi_FUN(thetahat.P, nuhat.P, Y, W)
A.P <- A_FUN(W)
A.Pinv <- solve(A.P)

rho.grad <- -t(t(xi) %*% A.Pinv %*% t(psi.P))
rho.FPT <- 1 * t(t(xi) %*% t(psi.P))

thresholds <- X[-length(X),] + diff(X)/2
thresholds <- tail(head(thresholds, -50),-50)
thresholds <- thresholds[seq(1, length(thresholds), by = 10)]

idx <- lapply(thresholds, function(thresh) X <= thresh)
C.list <- lapply(idx, function(i) list(C1 = which(i), C2 = which(!i)))

#--- Compute fitted versions of \theta
theta.fits <- lapply(C.list, function(split) {
  lapply(split, function(Cj) {
    thetahat <- thetahatC_FUN(Cj, Y, W)
    thetatilde.grad <- thetatildeC_FUN(Cj, thetahat.P, rho.grad)
    thetatilde.FPT  <- thetatildeC_FUN(Cj, thetahat.P, rho.FPT)
    thetatilde.FPT2 <- thetatildeC_FUN(Cj, thetahat.P, eta.new * rho.FPT)
    list(hat = thetahat, 
         tilde.grad = thetatilde.grad, 
         tilde.FPT  = thetatilde.FPT,
         tilde.FPT2 = thetatilde.FPT2)
  }) |> 
    purrr::list_transpose() |>
    lapply(simplify2array)
}) |>
  purrr::list_transpose() |>
  lapply(simplify2array)

#--- Extract fits
thetahat.C1 <- t(theta.fits$hat[,"C1",])
thetatilde.C1.grad <- t(theta.fits$tilde.grad[,"C1",])
thetatilde.C1.FPT  <- t(theta.fits$tilde.FPT[,"C1",])
thetatilde.C1.FPT2 <- t(theta.fits$tilde.FPT2[,"C1",])

thetahat.C2 <- t(theta.fits$hat[,"C2",])
thetatilde.C2.grad <- t(theta.fits$tilde.grad[,"C2",])
thetatilde.C2.FPT  <- t(theta.fits$tilde.FPT[,"C2",])
thetatilde.C2.FPT2 <- t(theta.fits$tilde.FPT2[,"C2",])

#--- Compute criteria and maximizing indices
crit.coef <- sapply(C.list, function(splits) length(splits$C1) * length(splits$C1)/n^2)
Delta <- crit.coef * apply(thetahat.C1 - thetahat.C2, 1, function(x) sum(x^2))
Deltatilde.grad <- crit.coef * apply(thetatilde.C1.grad - thetatilde.C2.grad, 1, function(x) sum(x^2))
Deltatilde.FPT  <- crit.coef * apply(thetatilde.C1.FPT - thetatilde.C2.FPT, 1, function(x) sum(x^2))
Deltatilde.FPT2 <- crit.coef * apply(thetatilde.C1.FPT2 - thetatilde.C2.FPT2, 1, function(x) sum(x^2))

Delta.list <- list(
  target = Delta,
  grad   = Deltatilde.grad,
  FPT    = Deltatilde.FPT,
  FPT2   = Deltatilde.FPT2
)

#--------------------------------------------------
#----- Plots
#--------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

my.cols <- gg_color_hue(4)

df.Delta <- mapply(
  function(dlt, id) bind_cols(threshold = thresholds, Delta = dlt, id = id),
  Delta.list, 
  list("target", "grad", "FPT", "FPT2"), 
  SIMPLIFY = FALSE) |>
  bind_rows() |>
  mutate(id = factor(id, levels = c("target", "grad", "FPT", "FPT2")))

Delta.max <- sapply(Delta.list, function(dlt) {
  c("x" = thresholds[which.max(dlt)], "y" = max(dlt))
})

plt <- ggplot(df.Delta, aes(x = threshold, y = Delta, color = id)) + 
  labs(y = expression("Criterion value"~Delta),
       x = expression("Threshold "*italic(t))) +
  ggtitle(expression("Criterion values"~Delta*"("*phantom("")*italic(C)[1]*","*italic(C)[2]*")"~"over candidate splits"~"{"*phantom("")*italic(C)[1]*","*italic(C)[2]*"}"),
          subtitle = expression(
            "\t"~italic(C)[1]~"= {"*phantom("")*italic(X)[italic(i)]*" : "*italic(X)[italic(i)] <= italic(t)*"}"~"and"~italic(C)[2]~"= {"*phantom("")*italic(X)[italic(i)]*" : "*italic(X)[italic(i)]*" > "*italic(t)*"}")
  ) + 
  geom_vline(xintercept = 0, linewidth = 1, color = "gray85") + 
  geom_hline(yintercept = 0, linewidth = 1, color = "gray85") + 
  geom_vline(xintercept = thresholds[which.max(Delta)], linewidth = 0.5, color = my.cols[1]) + 
  geom_vline(xintercept = thresholds[which.max(Deltatilde.grad)], linewidth = 0.5, color = my.cols[2]) +
  geom_vline(xintercept = thresholds[which.max(Deltatilde.FPT)], linewidth = 0.5, color = my.cols[3]) + 
  geom_vline(xintercept = thresholds[which.max(Deltatilde.FPT2)], linewidth = 0.5, color = my.cols[4]) + 
  geom_line(linewidth = 0.9) + 
  annotate("point", x = Delta.max["x", "target"], y = Delta.max["y", "target"], color = ggplot2::alpha(my.cols[1], 0.75), size = 2) + 
  annotate("point", x = Delta.max["x", "grad"], y = Delta.max["y", "grad"], color = ggplot2::alpha(my.cols[2], 0.75), size = 2) + 
  annotate("point", x = Delta.max["x", "FPT"], y = Delta.max["y", "FPT"], color = ggplot2::alpha(my.cols[3], 0.75), size = 2) + 
  annotate("point", x = Delta.max["x", "FPT2"], y = Delta.max["y", "FPT2"], color = ggplot2::alpha(my.cols[4], 0.75), size = 2) + 
  scale_color_manual(
    values = my.cols,
    labels = expression(Delta, 
                        tilde(Delta)^{grad}, 
                        tilde(Delta)[1]^{FPT}~eta~"="~1, 
                        tilde(Delta)[2]^{FPT}~eta~"="~sqrt(3/2))) + 
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.185,0.75)
  )
pdf("figures/example-criteria-comparison.pdf", width = 5.5, height = 4)
plt
dev.off()
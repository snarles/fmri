## approximating chi-squared sums sum beta_i X_i^2

library(pracma)

ss_samp <- function(betas, n) {
  p <- length(betas)
  colSums(betas * randn(p, n)^2)
}

kfunc0 <- function(betas, a) -1/(2 * length(betas)) * sum(1 - 2*a*betas)
kfunc1 <- function(betas, a) 1/(2 * length(betas)) * sum(1/(a - (1/(2 * betas))))
kfunc2 <- function(betas, a) -1/(2 * length(betas)) * sum((1/(a - (1/(2 * betas))))^2)

betas <- c(2, 1, 1)
ys <- ss_samp(betas, 1e5)
mean(ys)
var(ys)
2 * sum(betas^2)

a0 <- 0.1
delta <- -1e-4
kfunc0(betas, a0)
kfunc0(betas, a0 + delta)
kfunc0(betas, a0) + delta * kfunc1(betas, a0)
kfunc0(betas, a0) + delta * kfunc1(betas, a0) + (delta^2/2) * kfunc2(betas, a0)

kfunc1(betas, a0)
kfunc1(betas, a0 + delta)
kfunc1(betas, a0) + delta * kfunc2(betas, a0)



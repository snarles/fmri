## Generate data and try to guess mu, sigma2

library(pracma)

objective_function <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  sum((accs - rowMeans(repmat(t(pnorm(ts)), length(ks), 1)^repmat(t(t(ks-1)), 1, length(ts))))^2)
}

objective_grad <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^(ks-1))
  dmu <- rowMeans(pts^(ks-1) * (Ts - mu)/tau)
  dtau <- rowMeans(pts^(ks-1) * ((mu-Ts)^2 - tau)/(2 * tau^2))
  
}


## data settings
(mu <- 5 * rexp(1))
(sigma2 <- exp(rnorm(1)/2))
tau <- sigma2
n.reps <- 50
ks <- 2:100
K <- max(ks)

z_stars <- rnorm(n.reps, mu, sqrt(sigma2))
zmat <- pracma::randn(n.reps, K)
zmaxs <- t(apply(zmat, 1, cummax))[, ks - 1]
dim(zmaxs)
accs <- colMeans(zmaxs < z_stars)
plot(ks, accs, type = "l")

objective_function(accs, ks, mu, sigma2)
objective_function(accs, ks, mu - 0.1, sigma2 + 0.1)

numDeriv::grad(function(x) objective_function(accs, ks, x[1], x[2]),
               c(mu, sigma2))
numDeriv::hessian(function(x) objective_function(accs, ks, x[1], x[2]),
               c(mu, sigma2))
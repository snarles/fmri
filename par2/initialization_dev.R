## heuristics for starting values?

source("par2/objective_function.R")

generate_accs <- function(mu, tau, ks, n.reps = 100) {
  K <- max(ks)
  z_stars <- rnorm(n.reps, mu, sqrt(tau))
  zmat <- pracma::randn(n.reps, K)
  zmaxs <- t(apply(zmat, 1, cummax))[, ks - 1]
  dim(zmaxs)
  accs <- colMeans(zmaxs < z_stars)
  accs
}

## data settings
(mu <- 5 * rexp(1))
(tau <- exp(rnorm(1)/2))
n.reps <- 100
ks <- 2:100
accs <- generate_accs(mu, tau, ks, n.reps)
plot(ks, accs, type = "l")

## initialization
(ks.crit <- ks[which.min((accs - 0.5)^2)])
mu.init <- qnorm(1 - 1/(ks.crit - 1))
tau.init <- 1 ## just start with 1

res <- nlm(function(x) objective_function(accs, ks, x[1], x[2]),
    c(mu.init, tau.init))
(mu.hat <- res$estimate[1])
(tau.hat <- res$estimate[2])

list(of_orig = objective_function(accs, ks, mu, tau),
     of_opt = objective_function(accs, ks, mu.hat, tau.hat))

list(mu.init = mu.init, mu.hat = mu.hat)
list(tau.init = tau.init, tau.hat = tau.hat)

# Compute Pr[N(mu, tau) < max_{k-1} Z]

library(numDeriv)

acc_k <- function(k, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  mean(pnorm(ts)^(k-1))
}

acc_k_mu <- function(k, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  mean(pnorm(ts)^(k-1) * (ts - mu)/tau)
}

acc_k_tau <- function(k, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  mean(pnorm(ts)^(k-1) * ((mu-ts)^2 - tau)/(2 * tau^2))
}


help("numericDeriv")

k = 5; mu = rnorm(1); tau = rexp(1);
acc_k_mu(k, mu, tau)

grad(function(mu) acc_k(k, mu, tau), mu)
acc_k_mu(k, mu, tau)

grad(function(tau) acc_k(k, mu, tau), tau)
acc_k_tau(k, mu, tau)

hessian(function(x) acc_k(k, x[1], x[2]), c(mu, tau))

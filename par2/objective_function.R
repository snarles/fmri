## Objective function for squared-error minimization of 2-par curve

require(pracma)

##function itself
par2_acc_k <- function(ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  fs
}

##extrapolation
par2_extrapolate <- function(ks, accs, Ktarg, mc.reps = 1e4) {
  mu.init <- init.mu(accs, ks)
  tau.init <- 1
  res <- nlm(function(x) objective_function(accs, ks, x[1], x[2], mc.reps),
             c(mu.init, tau.init))
  (mu.hat <- res$estimate[1])
  (tau.hat <- res$estimate[2])
  par2_acc_k(Ktarg, mu.hat, tau.hat, mc.reps)
}

## objectives
par2_objective_function <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  sum((accs - fs)^2)
}

par2_objective_grad <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  dmu <- rowMeans(pts^Ks * (Ts - mu)/tau)
  dtau <- rowMeans(pts^Ks * ((mu-Ts)^2 - tau)/(2 * tau^2))
  c(sum(2 * (fs - accs) * dmu),
    sum(2 * (fs - accs) * dtau))
}

par2_objective_hess <- function(accs, ks, mu, tau, mc.reps = 1e4) {
  xs <- ((1:mc.reps) - 0.5)/mc.reps
  ts <- qnorm(xs, mean = mu, sd = sqrt(tau))
  Ts <- repmat(t(ts), length(ks), 1)
  pts <- repmat(t(pnorm(ts)), length(ks), 1)
  Ks <- repmat(t(t(ks-1)), 1, length(ts))
  fs <- rowMeans(pts^Ks)
  dmu <- rowMeans(pts^Ks * (Ts - mu)/tau)
  dtau <- rowMeans(pts^Ks * ((mu-Ts)^2 - tau)/(2 * tau^2))
  d2mu <- rowMeans(pts^Ks * (-1/tau + ((Ts - mu)/tau)^2))
  d2mu_tau <- rowMeans(pts^Ks *((mu - Ts)/(tau^2) + (Ts - mu)/tau * ((mu-Ts)^2 - tau)/(2 * tau^2)))
  d2tau <- rowMeans(pts^Ks * ((tau - 2 * (mu - Ts)^2)/(2 * tau^3) + (((mu-Ts)^2 - tau)/(2 * tau^2))^2))
  h11 <- 2 * sum(dmu^2 + (fs - accs) * d2mu)
  h12 <- 2 * sum(dmu*dtau + (fs - accs) * d2mu_tau)
  h22 <- 2 * sum(dtau^2 + (fs - accs) * d2tau)
  matrix(c(h11, h12, h12, h22), 2, 2)
}

attr(par2_objective_function, "gradient") <- par2_objective_grad
attr(par2_objective_function, "hessian") <- par2_objective_hess

## rough initialization
init.mu <- function(accs, ks) {
  (ks.crit <- ks[which.min((accs - 0.5)^2)])
  mu.init <- qnorm(1 - 1/(ks.crit - 1))
  mu.init
}
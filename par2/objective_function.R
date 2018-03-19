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
  if (min(accs) == 1) {
    return(rep(1, length(Ktarg)))
  }
  mu.init <- par2_initialize_mu(accs, ks)
  tau.init <- 1
  res <- tryCatch({
nlm(function(x) par2_objective_function(accs, ks, x[1], x[2]),
        c(mu.init, tau.init))},
    error = function(e) e)
  if ("error" %in% class(res)) {
    lala <- tryCatch({
      nlm(function(x) par2_objective_function(accs, ks, x[1], x[2]),
          c(mu.init, tau.init))},
      warning = function(e) e)
    if("sqrt(tau)" %in% lala) {
      res2 <- optimize(function(x) par2_objective_function(accs, ks, x, 0), interval = c(-20, 20))
      (mu.hat <- res2$minimum)
      tau.hat <- 0      
    } else {
      print(lala)
      stop("unknown error")
    }    
  } else {
    (mu.hat <- res$estimate[1])
    (tau.hat <- res$estimate[2])
  }
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
subroutine_infer_mu_tau <- function(mus, taus, ps) {
  p2norm <- function(mu, tau, mu2, tau2) pnorm((mu - mu2)/sqrt(tau + tau2))
  (qs2 <- qnorm(ps)^2)
  aa <- (qs2[2] - qs2[1])
  bb <- -2*(qs2[2]*mus[1] - qs2[1]*mus[2])
  cc <- qs2[2]*mus[1]^2 - qs2[1]*mus[2]^2 - qs2[1]*qs2[2] * (taus[1]-taus[2])
  if (bb < 4 * aa * cc) {
    return(c(subroutine_infer_mu(mus[2], taus[2], ps[2]),1))
  }
  (sols.mu <- c((-bb + sqrt(bb^2 - 4*aa*cc))/(2 * aa), (-bb - sqrt(bb^2 - 4*aa*cc))/(2 * aa)))
  (sols.tau <- (sols.mu - mus[1])^2/qs2[1] - taus[1])
  (sols.mu <- sols.mu[sols.tau > 0])
  (sols.tau <- sols.tau[sols.tau > 0])
  if (length(sols.mu) == 0) {
    return(c(subroutine_infer_mu(mus[2], taus[2], ps[2]),1))
  }
  sol.check <- sapply(1:length(sols.mu), 
                      function(j) 
                        sum(abs(sapply(1:2, function(i) p2norm(sols.mu[j], sols.tau[j], mus[i], taus[i])) - ps)))
  (mu.hat <- sols.mu[which.min(sol.check)])
  (tau.hat <- sols.tau[which.min(sol.check)])
  c(mu.hat, tau.hat)  
}

subroutine_infer_mu <- function(mu, tau, pp, tau0 = 1) {
  q2 <- qnorm(pp)^2
  aa <- 1
  bb <- -2 * mu
  cc <- mu^2 - (tau + tau0) * q2
  (sols.mu <- c((-bb + sqrt(bb^2 - 4*aa*cc))/(2 * aa), (-bb - sqrt(bb^2 - 4*aa*cc))/(2 * aa)))
  p2norm <- function(mu, tau, mu2, tau2) pnorm((mu - mu2)/sqrt(tau + tau2))
  sols.check <- abs(sapply(sols.mu, function(m) p2norm(m, tau0, mu, tau)) - pp)
  sols.mu[which.min(sols.check)]
}

approx_gaussian_maxk_moments <- function(k) {
  if (k == 1) return(c(0, 1))
  mu.coefs <- c(-0.676, 1.05)
  var.coefs <- c(0.048984, 1.15798, 0.001205)
  var.coefs_sm <- c(0.990, 0.648, 0.078)
  mu.hat <- mu.coefs[1] + mu.coefs[2] * sqrt(2 * log(k))
  if (k < 100) {
    pres.hat <- var.coefs_sm[1] + var.coefs_sm[2] * log(k) + var.coefs_sm[3] * log(k)^2
  } else {
    pres.hat <- var.coefs[1] + var.coefs[2] * log(k) + var.coefs[3] * log(k)^2
  }
  c(mu.hat, 1/pres.hat)
}

par2_initialize_mu <- function(accs, ks) {
  (ind <- which.min((accs - 0.5)^2))
  k <- ks[ind]
  mm <- approx_gaussian_maxk_moments(k-1)
  subroutine_infer_mu(mm[1], mm[2], accs[ind])
}
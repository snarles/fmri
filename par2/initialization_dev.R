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
(mu <- 3 * rexp(1))
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


####
##  Suppose
##   Pr[N(mu, tau) > N(mu_i, tau_i)] = p_i
##  for i = 1,2.  Find (mu, tau)
####

(mus <- rnorm(2))
(taus <- rexp(2))

(mu0 <- rnorm(1))
(tau0 <- rexp(1))

p2norm <- function(mu, tau, mu2, tau2) pnorm((mu - mu2)/sqrt(tau + tau2))
ps <- sapply(1:2, function(i) p2norm(mu0, tau0, mus[i], taus[i]))
ps

(qs2 <- qnorm(ps)^2)

i <- 1; (mu0 - mus[i])^2/qs2[i] - taus[i]
i <- 2; (mu0 - mus[i])^2/qs2[i] - taus[i]


(qs2[2] - qs2[1]) * mu0^2 - 
  2*(qs2[2]*mus[1] - qs2[1]*mus[2]) * mu0 + 
  qs2[2]*mus[1]^2 - qs2[1]*mus[2]^2 - qs2[1]*qs2[2] * (taus[1]-taus[2])

aa <- (qs2[2] - qs2[1])
bb <- -2*(qs2[2]*mus[1] - qs2[1]*mus[2])
cc <- qs2[2]*mus[1]^2 - qs2[1]*mus[2]^2 - qs2[1]*qs2[2] * (taus[1]-taus[2])

mu0
(sols.mu <- c((-bb + sqrt(bb^2 - 4*aa*cc))/(2 * aa), (-bb - sqrt(bb^2 - 4*aa*cc))/(2 * aa)))
(sols.tau <- (sols.mu - mus[1])^2/qs2[1] - taus[1])
(sols.mu <- sols.mu[sols.tau > 0])
(sols.tau <- sols.tau[sols.tau > 0])

sol.check <- sapply(1:length(sols.mu), 
                   function(j) 
                     sum(abs(sapply(1:2, function(i) p2norm(sols.mu[j], sols.tau[j], mus[i], taus[i])) - ps)))
sol.check


(mu.hat <- sols.mu[which.min(sol.check)])
(tau.hat <- sols.tau[which.min(sol.check)])

c(mu.hat, mu0)
c(tau.hat, tau0)


subroutine_infer_mu_tau <- function(mus, taus, ps) {
  p2norm <- function(mu, tau, mu2, tau2) pnorm((mu - mu2)/sqrt(tau + tau2))
  (qs2 <- qnorm(ps)^2)
  aa <- (qs2[2] - qs2[1])
  bb <- -2*(qs2[2]*mus[1] - qs2[1]*mus[2])
  cc <- qs2[2]*mus[1]^2 - qs2[1]*mus[2]^2 - qs2[1]*qs2[2] * (taus[1]-taus[2])
  (sols.mu <- c((-bb + sqrt(bb^2 - 4*aa*cc))/(2 * aa), (-bb - sqrt(bb^2 - 4*aa*cc))/(2 * aa)))
  (sols.tau <- (sols.mu - mus[1])^2/qs2[1] - taus[1])
  (sols.mu <- sols.mu[sols.tau > 0])
  (sols.tau <- sols.tau[sols.tau > 0])
  
  sol.check <- sapply(1:length(sols.mu), 
                      function(j) 
                        sum(abs(sapply(1:2, function(i) p2norm(sols.mu[j], sols.tau[j], mus[i], taus[i])) - ps)))
  (mu.hat <- sols.mu[which.min(sol.check)])
  (tau.hat <- sols.tau[which.min(sol.check)])
  c(mu.hat, tau.hat)  
}


(mus <- rnorm(2))
(taus <- rexp(2))

(mu0 <- rnorm(1))
(tau0 <- rexp(1))

p2norm <- function(mu, tau, mu2, tau2) pnorm((mu - mu2)/sqrt(tau + tau2))
ps <- sapply(1:2, function(i) p2norm(mu0, tau0, mus[i], taus[i]))
subroutine_infer_mu_tau(mus, taus, ps)

## how well does it work when using max of normals ??

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

(mu0 <- 2 + rnorm(1))
(tau0 <- rexp(1))
ks <- c(3, 4)
(ps <- par2_acc_k(ks, mu0, tau0))
(lala <- sapply(ks - 1, approx_gaussian_maxk_moments))
subroutine_infer_mu_tau(lala[1, ], lala[2, ], ps)

## assume tau = 1

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

(mu0 <- 2 + rnorm(1))
(tau0 <- 1)
k <- 4
(pp <- par2_acc_k(k, mu0, tau0))
(mms <- approx_gaussian_maxk_moments(k - 1))
subroutine_infer_mu(mms[1], mms[2], pp)

## function with backup

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

(mu0 <- 4 + rnorm(1))
(tau0 <- 0.5 + rexp(1))
ks <- c(30, 440)
(ps <- par2_acc_k(ks, mu0, tau0))
(lala <- sapply(ks - 1, approx_gaussian_maxk_moments))
subroutine_infer_mu_tau(lala[1, ], lala[2, ], ps)

####
##  back to initialization dev
####

## data settings
(mu <- 3 * rexp(1))
(tau <- exp(rnorm(1)/2))
n.reps <- 100
ks <- 2:100
accs <- generate_accs(mu, tau, ks, n.reps)
plot(ks, accs, type = "l")

## choose indices
(ind2 <- which.min((accs - 0.5)^2))
(ind1 <- which.min((accs - (1 + accs[ind2])/2)^2))
#stopifnot(ind1 < ind2)
inds <- c(ind1, ind2)
kz <- ks[inds]
ps <- accs[inds]
(lala <- sapply(kz - 1, approx_gaussian_maxk_moments))
#pars <- subroutine_infer_mu_tau(lala[1, ], lala[2, ], ps)
#pars
subroutine_infer_mu(lala[1, 2], lala[2, 2], ps[2])



(ind <- which.min((accs - 0.5)^2))
k <- ks[ind]
mm <- approx_gaussian_maxk_moments(k-1)
subroutine_infer_mu(mm[1], mm[2], accs[ind])

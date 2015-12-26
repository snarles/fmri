library(AlgDesign)
source("info_theory_sims/sim3source.R")

###
##  doing laplace approximation involve regularized logistic regression!
##  y ~ logistic(t(B) %*% x), st ||x||_2 small
###

str2vec <- function(s) {
  as.numeric(substring(s, seq(1,nchar(s),1), seq(1,nchar(s),1)))
}

eta <- function(x) log(1 + exp(x))
deta <- function(x) 1/(1 + exp(-x))
d2eta <- function(x) exp(-x)/(1 + exp(-x))^2

eta(1)
numDeriv::grad(eta, 1)
deta(1)
numDeriv::hessian(eta, 1)
d2eta(1)


pr_grid <- function(y, Bmat, reso = 5) {
  Btilde <- Bmat %*% diag(-y)
  xs <- as.matrix(gen.factorial(rep(reso, p))/reso * 2 * sqrt(2 * log(reso)))
  delta <- xs[2, 1] - xs[1, 1]
  mus <- xs %*% Btilde
  nms <- rowSums(xs^2)/2+ rowSums(eta(mus))
  delta^p * (1/sqrt(2*pi))^p * dim(xs)[1] * meanexp(-nms)
}

pr_laplace <- function(y, Bmat, log.p = FALSE) {
  p <- dim(Bmat)[1]; q <- dim(Bmat)[2]
  Bt <- Bmat %*% diag(-y)
  x0 <- opt_nll(y, Bmat)
  mu <- as.numeric(t(Bt) %*% x0)
  dn <- deta(mu); d2n <- d2eta(mu)
  (l0 <- nll(y, x0, Bmat))
#   gd <- x0 + Bt %*% dn
  hs <- eye(p) + Bmat %*% diag(d2n) %*% t(Bmat)
#   as.numeric((1/sqrt(2*pi))^p * 
#                exp(-l0 + 1/2 * t(gd) %*% solve(hs, gd)) * 
#                sqrt(det(2 * pi * solve(hs))))
  if (log.p) return(-l0 - log(det(hs))/2)
  exp(-l0)/sqrt(det(hs))
}


nll <- function(y, x, Bmat) {
  sum(x^2)/2 + sum(eta(-y * (t(Bmat) %*% x)))
}

opt_nll <- function(y, Bmat) {
  p <- dim(Bmat)[1]
  res <- optim(rep(0, p), function(x) nll(y, x, Bmat))
  res$par
}

prox_nll <- function(y, x, Bmat, x0 = opt_nll(y, Bmat)) {
  delta <- x - x0
  mu <- as.numeric(-y * (t(Bmat) %*% x0))
  dn <- deta(mu); d2n <- d2eta(mu)
  nll(y, x0, Bmat) + sum(delta * x0) + sum(delta^2)/2 +
    sum(dn * (-y * (t(Bmat) %*% delta))) + sum(d2n * (-y * (t(Bmat) %*% delta))^2)/2
}

p <- 2
q <- 3
Bmat <- 200 * randn(p, q)
mc.reps <- 1
X <- randn(mc.reps, p)
ps <- 1/(1 + exp(-X %*% Bmat))
Y <- (rand(mc.reps, q) < ps) + 0
y <- 2 * as.numeric(Y) - 1
pr_true <- pr_grid(y, Bmat, 300)
pr_true2 <- pr_grid(y, Bmat, 400)
pr_the <- pr_laplace(y, Bmat)
c(pr_true = pr_true, pr_true2 = pr_true2, pr_the = pr_the)

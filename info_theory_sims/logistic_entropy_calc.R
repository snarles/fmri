library(AlgDesign)
source("info_theory_sims/sim3source.R")

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

pr_laplace <- function(y, Bmat) {
  p <- dim(Bmat)[1]; q <- dim(Bmat)[2]
  Bt <- Bmat %*% diag(-y)
  x0 <- opt_nll(y, Bmat)
  mu <- as.numeric(t(Bt) %*% x0)
  dn <- deta(mu); d2n <- d2eta(mu)
  (l0 <- nll(y, x0, Bmat))
#   gd <- x0 + Bt %*% dn
  hs <- eye(p) + Bt %*% diag(d2n) %*% t(Bt)
#   as.numeric((1/sqrt(2*pi))^p * 
#                exp(-l0 + 1/2 * t(gd) %*% solve(hs, gd)) * 
#                sqrt(det(2 * pi * solve(hs))))
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
Bmat <- randn(p, q)
mc.reps <- 1e6
X <- randn(mc.reps, p)
ps <- 1/(1 + exp(-X %*% Bmat))
Y <- (rand(mc.reps, q) < ps) + 0
ylabs <- apply(Y, 1, function(v) paste(v, collapse = ""))
tab <- table(ylabs)
tab <- tab/sum(tab)

#(s <- "111")
(s <- sample(ylabs, 1))
y <- 2 * str2vec(s) - 1
(pr_emp <- tab[s])
(pr_true <- pr_grid(y, Bmat, 300))
(pr_the <- pr_laplace(y, Bmat))
# 
# x0 <- opt_nll(y, Bmat)
# x <- x0 + rnorm(p)/10
# nll(y, x, Bmat)
# prox_nll(y, x, Bmat, x0)
# 
# x0
# ds <- seq(-1, 1, 0.1)
# tvs <- sapply(ds, function(dd) nll(y, x + c(dd, 0), Bmat))
# pvs <- sapply(ds, function(dd) prox_nll(y, x + c(dd, 0), Bmat, x0))
# plot(ds, tvs, type = "l")
# lines(ds, pvs, col = "red")

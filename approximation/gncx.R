####
##  General non-central chi squared
####

library(pracma); library(magrittr)

rgchisq0 <- function(n, Sigma, mu) {
  p <- dim(Sigma)[1]
  L <- chol(Sigma)
  x <- L %*% (randn(p, n) + mu)
  colSums(x^2)
}

rgchisq <- function(n, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  x <- sqrt(ls) * (randn(p, n) + nu)
  colSums(x^2)
}

mgf_gchisq <- function(tt, Sigma, mu, log = FALSE) {
  p <- dim(Sigma)[1]
  tt <- t(t(tt))
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  tM <- repmat(tt, 1, p)
  lM <- repmat(t(ls), length(tt), 1)
  nM <- repmat(t(nu), length(tt), 1)
  if (log) {
    temp <- -1/2 * log(1 - 2*tM*lM) + (nM^2*lM*tM/(1 - 2*tM*lM))
    return(rowSums(temp))
  } else {
    temp <- 1/sqrt(1 - 2*tM*lM)*exp(nM^2*lM*tM/(1 - 2*tM*lM))
    return(apply(temp, 1, prod))
  }
}

## fourier inversion
finv <- function(fs, cs, x) {
  diffs <- Im(fs)[-1] - Im(fs)[-length(fs)]
  fM <- repmat(t(fs), length(x), 1)
  cM <- repmat(t(cs), length(x), 1)
  xM <- repmat(t(t(x)), 1, length(fs))
  vals <- exp(-fM * xM) * cM
  vals2 <- (vals[, -1, drop = FALSE] 
            + vals[, -length(fs), drop = FALSE])/2
  Re(rowSums(vals2 * diffs))/(2 * pi)
}

####
##  Upper bound using Markov
####

## exponential bound on Pr[<x]
mb_gchisq <- function(tt, x, Sigma, mu) {
  ms <- exp(mgf_gchisq(tt, Sigma, mu, TRUE))
  ms/exp(tt * x)
}

## derivative of log. expo. bound wrt t
prox_dleb <- function(tt, x, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  u <- 1/tt
  # actual
  sum(ls*(1 + nu^2)/(1-2*ls*tt)) + sum(2*nu^2*ls^2*tt/(1 - 2*ls*tt)^2) - x
  # second-order
  (-p/2)*u + sum((1+2*nu+3*nu^2)/ls)/4*u^2 - x
  # first-order
  (-p/2) * u - x
}

## exponential bound on log Pr[<x] plugging in t = -p/2x
# lmb_gchisq <- function(x, Sigma, mu) {
#   log(mb_gchisq(-.5*dim(Sigma)[1]/x, x, Sigma, mu))
# }

lmb_gchisq <- function(x, Sigma, mu) {
    p <- dim(Sigma)[1]
    res <- eigen(Sigma); ls <- res$values
    nu <- as.numeric(t(res$vectors) %*% mu)
    -1/2 * sum(log(1 + ls*p/x)) + p/2 * 
      (1 - sum(nu^2*ls/(1 + ls*p/x))/x)
}

qlmb_gchisq <- function(lprob, Sigma, mu, intv = c(1e-10, 1e3)) {
  ff <- function(x) (lmb_gchisq(x, Sigma, mu) - lprob)^2
  res <- optimise(ff, interval=intv)
  res$minimum
}

####
##  Lower bound using cap
####

log_vsph <- function(d, r=1)
  log(2/d) + (d/2)*log(pi) - lgamma(d/2) + d*log(r)
#c(exp(log_vsph(3)), 4/3 * pi)

## bound on log Pr[<x], creates function handle
cap_lb_ <- function(x, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  ## Compute height constants, etc
  nun <- nu/sqrt(f2(nu))
  nSn <- sum(nun^2/ls)
  #c(nSn, t(mu) %*% solve(Sigma, mu)/f2(mu))
  zh <- sqrt(x/nSn) * (nun/ls)
  h <- sqrt(x * nSn)
  m0 <- sqrt(f2(nu, zh))
  ## Volume of ellipse
  vc <- log_vsph(p) + (p/2)*log(x) - (1/2)*sum(log(ls))
  #c(h, sum(zh * nun))
  ## Function handle
  ff <- function(u) {
    ## volume of cap
    cv <- vc + pbeta(u, (p-1)/2, (p-1)/2, log.p=TRUE)
    ml <- m0 + sqrt(x)/sqrt(min(ls)) * sqrt(1-(1-u)^2) + 2*u*h
    dens <- -(p/2) * log(2*pi) - (ml^2)/2
    dens + cv
  }
  ff
}



####
##  Tests
####

# p <- 2;
# mu <- rnorm(p)
# Sigma <- cov(randn(2*p, p))
# #s1 <- rgchisq0(1e6, Sigma, mu)
# s1 <- rgchisq(1e6, Sigma, mu)
# cdf <- function(x) sum(s1 < x)/length(s1)
# cap_par <- 0.02
# s1[50] %>% {c(cap_lb_(., Sigma, mu)(cap_par), log(cdf(.)), lmb_gchisq(., Sigma, mu))}
# 
# x <- s1[50]
# ff <- cap_lb_(x, Sigma, mu)
# (1:100/1000) %>% plot(., sapply(., ff), type = "l")

# #s2 <- rgchisq(1e6, Sigma, mu)
# #cdf0 <- function(x) mean(s1 < x)
# #plot(sort(s1), sort(s2), type = 'l')
# #abline(0, 1, col = "red")
# #.2 %>%{c(mean(exp(.*s1)),mean(exp(.* s2)),mgf_gchisq(.,Sigma,mu))}
# fs <- seq(-3, 3, by = 0.002) * 1i
# cs <- mgf_gchisq(fs, Sigma, mu)
# #plot(Im(fs), Re(cs), type = "l")
# #plot(Im(fs), Im(cs), type = "l")
# xs <- seq(0, 1, length.out = 1e3)
# ps <- finv(fs, cs, xs)
# #10 %>% {c(sum(ps[xs < .])/sum(ps), mean(s1 < .))}
# cdf <- function(x) sum(ps[xs < x]) * (xs[2]-xs[1])
# x <- .5
# tt <- -seq(0,4/x,length.out=100); pus <- mb_gchisq(tt, x, Sigma, mu)
# plot(tt, log(pus), type = "l", ylim = c(log(cdf(x)) - 1, 0))
# abline(log(cdf(x)), 0)
# c(lmb_gchisq(x,Sigma,mu), log(min(pus)), log(cdf(x)))
# 
# ## test derivatives of log expo bound
# x <- 0.1
# tt0 <- -5
# numDeriv::grad(function(tt) mgf_gchisq(tt, Sigma, mu, TRUE) - tt*x, tt0)
# prox_dleb(tt0, x, Sigma, mu)

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
##  Tests
####

# mu <- rnorm(5)
# Sigma <- cov(randn(10, 5))
# #s1 <- rgchisq0(1e6, Sigma, mu)
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

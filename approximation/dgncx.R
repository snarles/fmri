####
##  DOUBLY General non-central chi squared
##  ie Z^T Omega Z, where Z ~ N(mu, Sigma)
####

library(pracma); library(magrittr)
library(MASS); library(lineId)

rgchisq <- function(n, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  x <- sqrt(ls) * (randn(p, n) + nu)
  colSums(x^2)
}

rgchisq0 <- function(n, Sigma, mu) {
  p <- dim(Sigma)[1]
  L <- chol(Sigma)
  x <- L %*% (randn(p, n) + mu)
  colSums(x^2)
}

rdgchisq <- function(n, Sigma, mu, Omega) {
  p <- dim(Sigma)[1]
  L <- chol(Omega)
  x <- mvrnorm(n, mu = mu, Sigma = Sigma)
  rowSums((x %*% L)^2)  
}

rdgchisq0 <- function(n, Sigma, mu, Omega) {
  si <- isqrtm(Sigma)
  sh <- sqrtm(Sigma)
  rgchisq0(n, sh %*% Omega %*% sh, (si %*% mu)[, 1])
}

####
##  Tests
####

# p <- 10; n<- 1e6
# Sigma <- cov(randn(2 * p, p))
# Omega <- cov(randn(2 * p, p))
# mu <- rnorm(p)
# s1 <- rdgchisq(n, Sigma, mu, Omega)
# s2 <- rdgchisq0(n, Sigma, mu, Omega)
# plot(sort(s1), sort(s2), type = "l")
# abline(0, 1, col ="red")

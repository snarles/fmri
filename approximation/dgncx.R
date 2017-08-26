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
  rowSums((x %*% t(L))^2)  
}

rdgchisq0 <- function(n, Sigma, mu, Omega) {
  si <- isqrtm(Sigma)
  sh <- sqrtm(Sigma)
  rgchisq0(n, sh %*% Omega %*% sh, (si %*% mu)[, 1])
}

TR <- function(A) sum(diag(A))
moments_dgchisq <- function(Sigma, mu, Omega) {
  mu2 <- mu %*% t(mu)
  m <- TR(Omega %*% (Sigma + mu2))
  m2 <- (t(mu) %*% Omega %*% mu)^2 + 
    4 * t(mu) %*% Omega %*% Sigma %*% Omega %*% mu +
    2 * TR(Omega %*% Sigma %*% Omega %*% Sigma) + 
    (TR(Omega %*% Sigma))^2 + 
    2 * (t(mu) %*% Omega %*% mu) * TR(Omega * Sigma)
  v <- 2 * TR(Omega %*% Sigma %*% Omega %*% (Sigma + 2 * mu2))
  list(m = m, m2 = m2[1], v = v[1])
}

####
##  Tests
####

p <- 5; n<- 1e6
Sigma <- cov(randn(2 * p, p))
Omega <- cov(randn(2 * p, p))
mu <- rnorm(p)
s1 <- rdgchisq(n, Sigma, mu, Omega)
s2 <- rdgchisq0(n, Sigma, mu, Omega)
plot(sort(s1), sort(s2), type = "l")
abline(0, 1, col ="red")

moments_dgchisq(Sigma, mu, Omega)
list(m = mean(s1), m2 = mean(s1^2), v = var(s1))
list(m = mean(s2), m2 = mean(s2^2), v = var(s2))

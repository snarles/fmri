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

####
##  Tests
####

mu <- rnorm(5)
Sigma <- cov(randn(10, 5))
s1 <- rgchisq0(1e6, Sigma, mu)
s2 <- rgchisq(1e6, Sigma, mu)
plot(sort(s1), sort(s2), type = 'l')
abline(0, 1, col = "red")

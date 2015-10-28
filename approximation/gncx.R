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

mgf_gchisq <- function(t, Sigma, mu, log = FALSE) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  ans <- -1/2 * sum(log(1 - 2*t*ls)) + sum(nu^2*ls*t/(1 - 2*t*ls))
  if (log) return(ans)
  exp(ans)
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

-1 %>%{c(mean(exp(.*s1)),mean(exp(.* s2)),mgf_gchisq(.,Sigma,mu))}

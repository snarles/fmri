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
##  Tests
####

mu <- rnorm(5)
Sigma <- cov(randn(10, 5))
s1 <- rgchisq0(1e6, Sigma, mu)
s2 <- rgchisq(1e6, Sigma, mu)
plot(sort(s1), sort(s2), type = 'l')
abline(0, 1, col = "red")

.2 %>%{c(mean(exp(.*s1)),mean(exp(.* s2)),mgf_gchisq(.,Sigma,mu))}

fs <- (seq(-2, 2, by = 0.005) + 0.005 * runif(801)) * 1i
cs <- mgf_gchisq(fs, Sigma, mu)
xs <- seq(0, 80, by = 0.5)
ps <- finv(fs, cs, xs)
10 %>% {c(sum(ps[xs < .])/sum(ps), mean(s1 < .))}

fs <- 5 * fs
cs <- 1/(1 - 2*fs)^(p/2)


plot(Im(fs), Re(cs), type = "l")
plot(Im(fs), Im(cs), type = "l")
finv(fs, cs, 1)
dchisq(1:10, p)



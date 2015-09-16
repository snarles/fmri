####
##  Error for random two-class gaussian classification
####

library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc

## most naive implementation of mc()
mc2c <- function(Sigma, mc.reps = 1000) {
  K <- 2
  p <- dim(Sigma)[1]
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- mvrnorm(n = K, mu = rep(0, p), Sigma = Sigma)
                  ys <- mus + randn(K, p)
                  lbls <- knn(mus, ys, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

## Pr[||X||^2 < ||X - D||^2] where D ~ N(0, 2Lambda)
mc2c_1 <- function(Sigma, mc.reps = 1000) {
  p <- dim(Sigma)[1]  
  lambdas <- eigen(Sigma)$values
  xs <- randn(p, mc.reps)
  ds <- sqrt(2 * lambdas) * randn(p, mc.reps)
  xn <- colSums(xs^2)
  zn <- colSums((xs + ds)^2)
  sum(xn > zn)/mc.reps  
}

## 1 - E[Phi(sqrt(S)/2)], where S = sum 4lambda_i^2 Z_i^2
## Uses monte carlo
mc2c_2 <- function(Sigma, mc.reps = 1e4) {
  p <- dim(Sigma)[1]  
  lambdas <- eigen(Sigma)$values
  as <- (2 * lambdas)
  zs <- randn(p, mc.reps)
  ss <- colSums(as * (zs^2))
  1 - mean(pnorm(sqrt(ss)/2))
}

## 1 - E[Phi(sqrt(S)/2)], where S = sum 4lambda_i^2 Z_i^2
## Uses Satterthwaithe approximation
mc2c_3 <- function(Sigma, mc.reps = 1e4) {
  p <- dim(Sigma)[1]  
  lambdas <- eigen(Sigma)$values
  as <- (2 * lambdas)
  df <- sum(as)^2/sum(as^2)
  mu <- sum(as)
  seq <- ((1:mc.reps) - 0.5)/mc.reps
  samp <- mu/df * qchisq(seq, df = df)
  1 - mean(pnorm(sqrt(samp)/2))
}

## compares sum a_i Z_i^2 with satter approx
satter_test <- function(as, mc.reps = 1e4) {
  p <- length(as)
  zs <- randn(p, mc.reps)
  ss <- colSums(as * (zs^2))
  #hist(ss)
  df <- sum(as)^2/sum(as^2)
  mu <- sum(as)
  ss2 <- mu/df * rchisq(mc.reps, df = df)
  plot(sort(ss), sort(ss2));abline(a = 0, b=1, col = 'red')
}
satter_test(as, 1e5)

Sigma <- 0.1 * rchisq(1, 5) * cov(randn(40, 10))
mc2c(Sigma, 1e4)
mc2c_1(Sigma, 1e4)
mc2c_2(Sigma, 1e4)
mc2c_3(Sigma, 1e4)

p <- length(as)
lambdas <- eigen(Sigma)$values
as <- (2 * lambdas)^2
zs <- randn(p, mc.reps)
ss <- colSums(as * (zs^2))
c(mean(ss), sum(as))
c(var(ss), 2 * sum(as^2))



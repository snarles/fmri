library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')

## "Ideally the curve mc(Sigma, K) only depends on a one-dimensional
##  function of Sigma, i.e. mc(Sigma, K) = g(f(Sigma), K)
##  where f(Sigma) is real-valued."

## most naive implementation of mc()
mc <- function(Sigma, K, mc.reps = 1000) {
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

## gets entire curve of mc, using hypergeometric resampling
mcK <- function(Sigma, Kmax, N = 200, mc.reps = 1000) {
  if (mc.reps == 1) {
    p <- dim(Sigma)[1]
    mus <- mvrnorm(n = N, mu = rep(0, p), Sigma = Sigma)
    ys <- mus + randn(N, p)
    dm <- fastPdist2(ys, mus)
    nbads <- sapply(1:N, function(i) sum(dm[i,] < dm[i,i]))
    Km <- repmat(1:Kmax, N, 1)
    Bm <- repmat(t(t(nbads)), 1, Kmax)
    temp <- N - Bm - Km + 1
    temp2 <- temp
    temp2[temp < 1] <- 1
    miscs <- 1 - exp(lgamma(N - Bm) + lgamma(N - Km + 1) -
                       lgamma(N) - lgamma(temp2))
    miscs[temp < 1] <- 1
    return(colMeans(miscs))
  }
  else {
    lala <- sapply(1:mc.reps, function(i) mcK(Sigma, Kmax, N, 1))
    return(rowMeans(lala))
  }
}

mi <- function(Sigma) {
  SigmaY <- Sigma + eye(dim(Sigma)[1])
  -1/2 * (log(det(Sigma)) - log(det(SigmaY)))
}


# mean of exp(v)
meanexp <- function(v) {
  vm <- max(v)
  mean(exp(v - vm)) * exp(vm)
}

## High-dimensional asymptotic misclassification curve for Sigma = cc/sqrt(d) * eye(d)
mcK_I <- function(cc, maxK, mc.reps = 1e4) {
  samp <- qnorm(((1:mc.reps) - 0.5)/mc.reps)
  temp <- log(1 - pnorm(samp - sqrt(cc)))
  tmat <- ((1:maxK) - 1) * repmat(temp, maxK, 1)
  1 - apply(tmat, 1, meanexp)
}

####
##  Learning curves for multi-class gaussian classification
####

library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc

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
mcK <- function(Sigma, Kmax, N = 1000, mc.reps = 1000) {
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


Kmax <- 20
N <- 200
mc.reps <- 5000

expment1 <- function(seed) {
  set.seed(seed)
  p <- 5
  Sigma <- cov(randn(2*p, p))
  curve <- mcK(Sigma, Kmax, N, mc.reps)
  mutual <- mi(Sigma)
  list(curve = curve, mutual = mutual, Sigma = Sigma)
}

library(lineId)
res <- lclapply(1:20, expment1, mc.cores = 7)
curves <- do.call(cbind, curve)
plot(res$curve[[1]])
zattach(res)
mutual <- unlist(mutual)
cols <- rainbow(length(mutual))[order(mutual)]
par(bg = "grey")
plot(1:Kmax, 1:Kmax/Kmax, pch = ".", col = "grey")
for (ii in 1:length(Sigma)) {
  lines(1:Kmax, curve[[ii]], col = cols[[ii]])
}


plot(mutual, curves[10, ])
plot(curves[5, ], curves[10, ])


par(bg = "white")






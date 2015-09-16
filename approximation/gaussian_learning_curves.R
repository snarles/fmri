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


expment <- function(seed) {
  set.seed(seed)
  p <- 5
  Sigma <- cov(randn(2*p, p))
  curve <- mcK(Sigma, 20, 100, 1000)
  mutual <- mi(Sigma)
  list(curve = curve, mutual = mutual, Sigma = Sigma)
}


N <- 1e4
Sigma <- cov(randn(5, 10))
mc(Sigma, 10, 1e5)
.5/sqrt(1e5)
mcK(Sigma, 20, 1e4, 2)[10]

## timing experiments

Sigma <- cov(randn(10, 10))


#Nlvs <- rep(c(1e2, 1e3, 1e4), each = 3)
Nlvs <- c(100, 200, 500, 1000, 1100, 1200, 1300, 1400)
times <- Nlvs * 0
Kmax <- 20
datas <- matrix(0, length(Nlvs), Kmax)
mc.reps <- 1
for (ii in 1:length(Nlvs)) {
  t1 <- proc.time()
  datas[ii, ] <- mcK(Sigma, Kmax, Nlvs[ii], 100)
  t2 <- proc.time() - t1
  times[ii] <- t2[3]
}

gt <- colMeans(datas[-(1:2), ])
gt[20]
plot(Nlvs, (datas[, 10] - gt[10])^2)
1/(datas[, 10] - gt[10])^2
times
plot(times, 1/(datas[, 10] - gt[10])^2, ylim = c(0, 2e7), xlim = c(0, 10))
plot(times, Nlvs)
## conclusion: choose N = 200













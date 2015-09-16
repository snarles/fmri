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
mc.reps <- 1000
.5/sqrt(mc.reps)

expment1 <- function(seed) {
  set.seed(seed)
  Sigma <- rchisq(1, 5) * cov(randn(8, 5))
  curve <- mcK(Sigma, Kmax, N, mc.reps)
  mutual <- mi(Sigma)
  list(curve = curve, mutual = mutual, Sigma = Sigma)
}

expment2 <- function(seed) {
  set.seed(seed)
  Sigma <- scalings[seed] * Sigma0
  curve <- mcK(Sigma, Kmax, N, mc.reps)
  mutual <- mi(Sigma)
  list(curve = curve, mutual = mutual, Sigma = Sigma)
}

Sigma0 <- cov(randn(8, 5))
scalings <- (5:15)/10

library(lineId)
res <- lclapply(1:100, expment1, mc.cores = 7)
#res <- lclapply(1:length(scalings), expment2, mc.cores = 7)
zattach(res)
curves <- do.call(cbind, curve)
mutual <- unlist(mutual)
cols <- rainbow(length(mutual))[order(mutual)]
par(bg = "grey")
plot(1:Kmax, 1:Kmax/Kmax, pch = ".", col = "grey")
for (ii in 1:length(Sigma)) {
  lines(1:Kmax, curve[[ii]], col = cols[[ii]])
}


trs <- sapply(Sigma, function(s) sum(diag(s)))
tr2s <- sapply(Sigma, function(s) sum(diag(s^2)))
res <- lm(log(curves[2, ]) ~ trs + tr2s + trs^2 + tr2s^2)
plot(log(curves[2, ]), res$fitted.values)

plot(trs, log(curves[2, ]))
plot(tr2s, log(curves[2, ]))

plot(mutual, log(curves[2, ]))
plot(curves[2, ], curves[10, ])
plot(curves[2, ], curves[20, ])
plot(curves[10, ], curves[20, ])


par(bg = "white")






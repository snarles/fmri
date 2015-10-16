library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')


## mu ~ N(0, I)
## y ~ N(mu*, sigma^2 I)

## most naive implementation of mc() for identity cov
mc_ident <- function(p, sigma2, K, mc.reps = 1000) {
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- randn(K, p)
                  ys <- mus + sqrt(sigma2) * randn(K, p)
                  lbls <- knn(mus, ys, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

## uses noncentral chi squared
mc_ident2 <- function(p, sigma2, K, mc.reps = 1000) {
  alpha <- sigma2/(1 + sigma2)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y2 <- (1 + sigma2) * rchisq(1, df = p)
                  d1 <- alpha * rchisq(1, df = p, ncp = alpha * y2)
                  ds <- rchisq(K - 1, df = p, ncp = y2)
                  (min(ds) < d1)
                })
  mean(mcs)
}

## actually generates a gaussian instead of a chi-squared
rchisq_g <- function(n, df, ncp = 0) {
  mu <- df + ncp
  vv <- 2 * df + 4 * ncp
  sqrt(vv) * rnorm(n) + mu
}

## uses gaussian--approximate!!
mc_ident3 <- function(p, sigma2, K, mc.reps = 1000) {
  alpha <- sigma2/(1 + sigma2)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y2 <- (1 + sigma2) * rchisq_g(1, df = p)
                  if (y2 > 0) {
                    d1 <- alpha * rchisq_g(1, df = p, ncp = alpha * y2)
                    ds <- rchisq_g(K - 1, df = p, ncp = y2)
                    return(min(ds) < d1)
                  } else {
                    return(1 - (1/K)) ## default in case of error
                  }
                })
  mean(mcs)
}

mean(rchisq(1e3, 3, 2.2))
mean(rchisq_g(1e3, 3, 2.2))
var(rchisq(1e3, 3, 2.2))
var(rchisq_g(1e3, 3, 2.2))

## Gaussian approximation is unusably bad

mc_ident(10, 3, 10, 1e4)
mc_ident2(10, 3, 10, 1e4)
mc_ident3(10, 3, 10, 1e4)

mc_ident(20, 1, 100, 1e4)
mc_ident2(20, 1, 100, 1e4)
mc_ident3(20, 1, 100, 1e4)

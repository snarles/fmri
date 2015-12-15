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
## same Evalue as mc_ident
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




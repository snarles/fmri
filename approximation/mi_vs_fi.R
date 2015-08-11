####
## MUTUAL INFORMATION VS FISHER INFORMATION IN THE GAUSSIAN CASE
## (X, Y) ~ N(0, Sigma)
## Mutual information I(X; Y)
## Fisher information J(x) given Y
####

library(pracma)
library(magrittr)
source('transfer/source.R')
pX <- 20
pY <- 20
p <- pX + pY
Sigma <- 1/p/2 * randn(p, 2 * p) %>% {. %*% t(.)}
SigmaX <- Sigma[1:pX, 1:pX]
SigmaY <- Sigma[-(1:pX), -(1:pX)]
SigmaXY <- Sigma[1:pX, -(1:pX)]
## Matrix giving mean of Y conditional on X
Amat <- t(solve(SigmaX, SigmaXY))
## Variance of Y conditional on X
SigmaYcX <- SigmaY - t(SigmaXY) %*% solve(SigmaX, SigmaXY)
Rmat <- isqrtm(SigmaX) %*% SigmaXY %*% isqrtm(SigmaY)

####
## MUTUAL INFORMATION
## Claim: Mutual information is -1/2 log det(I - R^T R)
####

## Usual formula
mi0 <- (1/2) * (log(det(SigmaX)) + log(det(SigmaY)) - log(det(Sigma)))
## New formula
s <- svd(Rmat)$d
mi1 <- -1/2 * sum(log(1-s^2))
## Comparison
print(c(mi0, mi1))

####
## FISHER INFORMATION
## Claim: Fisher information = inv(SigmaX) * (SigmaXY * inv(SigmaYcX) * SigmaYX) * inv(SigmaX)
####

## Formula
fi0 <- solve(SigmaX, SigmaXY) %*% solve(SigmaYcX, t(SigmaXY)) %*% solve(SigmaX)

## BY gives an unbiased estimator for X given Y
g <- SigmaXY %*% solve(SigmaYcX, t(SigmaXY))
B <- SigmaX %*% solve(g, SigmaXY) %*% solve(SigmaYcX)

## Generate multivariate normals
library(MASS)
XY <- mvrnorm(n = 100, mu=rep(0, p), Sigma=Sigma)
X <- XY[, 1:pX]
Y <- XY[, -(1:pX)]
Xhat <- Y %*% B

## TODO: Confirm that they are unbiased for X, etc.

####
##  FISHER INFORMATION AND MUTUAL INFORMATION
##  Claim: 1/2 log det(I - J) and MI converge if SigmaYcX converges to SigmaY
##  I.e. this occurs if error size increases
####

## Generate matrices
experiment1 <- function(pX, pY, X_df, noise_df, noise_size) {
  SigmaX <- (1/pX/X_df) * randn(pX, X_df * pX) %>% {. %*% t(.)}
  B <- 1/pX * randn(pX, pY)
  SigmaXY <- SigmaX %*% B
  SigmaYX <- t(SigmaXY)
  SigmaE <- (noise_size/pY/noise_df) * randn(pY, noise_df * pY) %>% {. %*% t(.)}
  SigmaY <- t(B) %*% SigmaX %*% B + SigmaE
  SigmaYcX <- SigmaY - t(SigmaXY) %*% solve(SigmaX, SigmaXY)
  Rmat <- isqrtm(SigmaX) %*% SigmaXY %*% isqrtm(SigmaY)
  s <- svd(Rmat)$d
  mi <- -1/2 * sum(log(1-s^2))
  fi <- solve(SigmaX, SigmaXY) %*% solve(SigmaYcX, t(SigmaXY)) %*% solve(SigmaX)
  v <- eigen(fi)$values
  fi2mi <- -1/2 * sum(log(1 - v))
  c(mi, fi2mi)
}

#X_df <- 2
#noise_df <- 2
#noise_size <- 2

experiment1(20, 20, 2, 3, 2)
experiment1(20, 20, 2, 3, 3)
experiment1(20, 20, 2, 3, 5)
experiment1(20, 20, 2, 3, 10)

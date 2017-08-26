library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')

TR <- function(a) sum(diag(a))

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

## uses conditioning on Y
mc2 <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
                  y_mus <- randn(p, K-1) + y
                  nms <- apply(y_mus, 2, nmS)
                  #nms <- colSums(Ha_S %*% y_mus)
                  min(nms) < t(y_mu) %*% Sigma %*% y_mu
                })
  mean(mcs)  
}



## APPROXIMATE uses spherical Y, mu
mc3a <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(rnorm(p))
                  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
                  y_mus <- randn(p, K-1) + y
                  nms <- apply(y_mus, 2, nmS)
                  #nms <- colSums(Ha_S %*% y_mus)
                  min(nms) < rchisq(1, d, 0)
                })
  mean(mcs)  
}

## APPROXIMATE separates mu* and mu_i
mc3b <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mus <- randn(p, K-1) + y
                  nms <- apply(y_mus, 2, nmS)
                  #nms <- colSums(Ha_S %*% y_mus)
                  min(nms) < t(y_mu) %*% Sigma %*% y_mu
                })
  mean(mcs)  
}

## APPROXIMATE uses approximation for n0
mc3c <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  Omega <- solve(Sigma)
  Amat <- eye(p) - solve(eye(p) + Omega)
  Ha_Y <- lineId::sqrtm(eye(p) + Omega)
  Ha_mu <- lineId::sqrtm(Amat)
  nmS <- function(v) t(v) %*% Sigma %*% v
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y <- as.numeric(Ha_Y %*% rnorm(p))
                  y_mus <- randn(p, K-1) + y
                  nms <- apply(y_mus, 2, nmS)
                  en0 <- as.numeric(t(y)%*%Amat%*%Sigma%*%Amat%*%y + TR(Amat%*%Sigma))
                  vn0 <- as.numeric((4*t(y)%*%Amat%*%Sigma%*%Amat%*%Sigma%*%Amat%*%y+
                                       2 * TR(Amat%*%Sigma%%Amat%*%Sigma)))
                  (ss <- (-(2 * en0/d) + sqrt((2*en0/d)^2 + 2 * vn0/d))/2)
                  (ll <- en0/ss - d)
                  if (ll < 0) {print(paste("ll=", ll)); y_err <<- y}
                  n0 <- ss * rchisq(1, d, ll)
                  min(nms) < n0
                })
  mean(mcs)  
}

## Tried to make faster version??
# mc2 <- function(Sigma, K, mc.reps = 1000) {
#   p <- dim(Sigma)[1]
#   Omega <- solve(Sigma)
#   Amat <- eye(p) - solve(eye(p) + Omega)
#   
#   y <- mvrnorm(mc.reps, mu=rep(0,p), Sigma=eye(p)+Omega)
#   y_mu <- mvrnorm(mc.reps, mu=rep(0,p), Sigma=Amat)
#   y_mu <- y_mu + (y %*% Amat)
#   ps <- 1-(1-runif(mc.reps))^(1/(K-1))
#   yn <- rowSums(y^2)
#   
#   mean(mcs) 
# }

####
## Multivariate Bayesian Regression: Computational Aspects
####

## Simulated Data

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
f2 <- function(x) sum(x^2)

n <- 30
nte <- 20
pX <- 60 # number of X-features
pY <- 70 # number of Y-responses
sigma_b <- 0.1 # size of each random coefficent
# noise covariance
Sigma_e <- 10 * (1/10/pY) * randn(10 * pY, pY) %>% {t(.) %*% .}
# Generate X in a way that columns are correlated
Xall <- randn(n + nte, 20) %*% abs(randn(20, pX)) + randn(n + nte, pX)
X <- Xall[1:n, ]
Xte <- Xall[-(1:n), ]
# Generate coefficents
B <- sigma_b * randn(pX, pY)
Y <- X %*% B + mvrnorm(n = n, mu = rep(0, pY), Sigma = Sigma_e)
Yte <- Xte %*% B + mvrnorm(n = nte, mu = rep(0, pY), Sigma = Sigma_e)


## A basic ridge estimate of B (ignores correlations)

Bvec <- as.numeric(B)
(err_0 <- Norm(Yte)^2)
lambda <- mean(diag(Sigma_e))/(sigma_b^2)
ridge_mu <- solve(t(X) %*% X + lambda * eye(pX), t(X) %*% Y)
ridge_Yte <- Xte %*% ridge_mu
(err_ridge <- f2(Yte - ridge_Yte))

## Compute posterior mean and covariance: the naive way

Omega_e <- solve(Sigma_e)
IX <- eye(pY) %x% X
Sigma_b <- sigma_b^2 * eye(pX)
Omega_b <- 1/sigma_b^2 * eye(pX)
Sigma_B <- sigma_b^2 * eye(pX * pY)
Omega_B <- 1/sigma_b^2 * eye(pX * pY)
Sigma_E <- Sigma_e %x% eye(n)
Omega_E <- Omega_e %x% eye(n)
yVec <- as.numeric(Y)
xtx <- t(X) %*% X
xt.o.y <- t(IX) %*% Omega_E %*% yVec

#post_cov <- solve(Omega_e %x% xtx + Omega_B)
tt <- proc.time()
post_mu <- solve(Omega_e %x% xtx + Omega_B, xt.o.y)
proc.time() -tt
Yte_post <- Xte %*% matrix(post_mu, pX, pY)
(err_post <- f2(Yte - Yte_post))

## Compute posterior using SVD

res <- svd(rbind(X, matrix(0, pX - n, pX)))
V <- res$v
D <- diag(res$d)
f2(xtx - V %*% D^2 %*% t(V))

DU <- diag(res$d) %*% t(res$u)
IV <- eye(pY) %x% V
IDU <- eye(pY) %x% DU
icov_post <- Omega_e %x% xtx + Omega_B
icov_svd <- IV %*% (Omega_e %x% D^2 + (1/sigma_b^2) * eye(pY * pX)) %*% t(IV)
#f2(icov_post - icov_svd)
#cov_svd <- IV %*% solve(Omega_e %x% D^2 + (1/sigma_b^2) * eye(pY * pX), t(IV))
#f2(cov_svd - post_cov)
tt <- proc.time()
mu_svd <- IV %*% solve(Omega_e %x% D^2 + (1/sigma_b^2) * eye(pY * pX), t(IV) %*% t(IX) %*% Omega_E %*% yVec)
proc.time() -tt
Yte_svd <- Xte %*% matrix(mu_svd, pX, pY)
(err_svd <- f2(Yte - Yte_svd))

#f2(post_cov_svd - post_cov)
#max(abs((post_cov_svd - post_cov)))

## Avoid computing the entire inverse
temp0 <- (eye(pY) %x% Xte) %*% IV

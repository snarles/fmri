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
Sigma_b <- 0.01 * diag(abs(rnorm(pY))) # size of each random coefficent
# noise covariance
Sigma_e <- 10 * (1/10/pY) * randn(10 * pY, pY) %>% {t(.) %*% .}
# Generate X in a way that columns are correlated
Xall <- randn(n + nte, 20) %*% abs(randn(20, pX)) + randn(n + nte, pX)
X <- Xall[1:n, ]
Xte <- Xall[-(1:n), ]
# Generate coefficents
B <- mvrnorm(n = pX, mu = rep(0, pY), Sigma = Sigma_b)
Y <- X %*% B + mvrnorm(n = n, mu = rep(0, pY), Sigma = Sigma_e)
Yte <- Xte %*% B + mvrnorm(n = nte, mu = rep(0, pY), Sigma = Sigma_e)


## A basic ridge estimate of B (ignores correlations)

Bvec <- as.numeric(B)
(err_0 <- Norm(Yte)^2)
lambda <- mean(diag(Sigma_e))/mean(diag(Sigma_b))
ridge_mu <- solve(t(X) %*% X + lambda * eye(pX), t(X) %*% Y)
ridge_Yte <- Xte %*% ridge_mu
(err_ridge <- f2(Yte - ridge_Yte))

## Compute posterior mean and covariance: the naive way

Omega_e <- solve(Sigma_e)
IX <- eye(pY) %x% X
Sigma_B <- Sigma_b %x% eye(pX)
Omega_B <- diag(1/diag(Sigma_b)) %x% eye(pX)
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


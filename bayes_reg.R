####
## Multivariate Bayesian Regression: Computational Aspects
####

## Simulated Data

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
f2 <- function(x) sum(x^2)

n <- 20
pX <- 5 # number of X-features
pY <- 7 # number of Y-responses
sigma_b <- 0.1 # size of each random coefficent
# noise covariance
Sigma_e <- 10 * (1/10/pY) * randn(10 * pY, pY) %>% {t(.) %*% .}
# Generate X in a way that columns are correlated
X <- randn(n, 20) %*% abs(randn(20, pX)) + randn(n, pX)
# Generate coefficents
B <- sigma_b * randn(pX, pY)
E <- mvrnorm(n = n, mu = rep(0, pY), Sigma = Sigma_e)
Y <- X %*% B + E

## A basic ridge estimate of B (ignores correlations)

Bvec <- as.numeric(B)
Norm(Bvec)^2
lambda <- mean(diag(Sigma_e))/(sigma_b^2)
ridge_mu <- as.numeric(solve(t(X) %*% X + lambda * eye(pX), t(X) %*% Y))
(err_ridge <- Norm(Bvec - ridge_mu)^2)

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

#f2(Omega_e %x% xtx - t(IX) %*% Omega_E %*% IX)
# check precision
#Norm(Omega_E %*% yVec - as.numeric(t(solve(Sigma_e, t(Y)))))

post_cov <- solve(Omega_e %x% xtx + Omega_B)
post_mu <- solve(Omega_e %x% xtx + Omega_B, xt.o.y)
# check precision
#Norm((Omega_e %x% xtx + Omega_B) %*% post_mu - xt.o.y)
(err_post <- Norm(post_mu - Bvec)^2)


plot(Bvec, ridge_mu)

plot(Bvec, post_mu)
plot(ridge_mu, post_mu)

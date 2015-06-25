####
## Multivariate Bayesian Regression: Computational Aspects
####

## Simulated Data

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
f2 <- function(x) sum(x^2)
solvediag <- function(D, b) 1/diag(D) * b
#lala <- rnorm(10); dada <-abs(diag(rnorm(10)))
#solve(dada, lala)
#solvediag(dada, lala)


n <- 5
nte <- 20
pX <- 10 # number of X-features
pY <- 10 # number of Y-responses
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
Omega_b <- diag(1/diag(Sigma_b))
Sigma_B <- Sigma_b %x% eye(pX)
Omega_B <- diag(1/diag(Sigma_b)) %x% eye(pX)
Sigma_E <- Sigma_e %x% eye(n)
Omega_E <- Omega_e %x% eye(n)
yVec <- as.numeric(Y)
xtx <- t(X) %*% X
xt.o.y <- t(IX) %*% Omega_E %*% yVec

tt <- proc.time()
icov <- Omega_e %x% xtx + Omega_B
post_cov <- solve(Omega_e %x% xtx + Omega_B)
proc.time() -tt
tt <- proc.time()
post_mu <- solve(Omega_e %x% xtx + Omega_B, xt.o.y)
proc.time() -tt
Yte_post <- Xte %*% matrix(post_mu, pX, pY)
(err_post <- f2(Yte - Yte_post))

## Using simultaneous diagonalization of Omega_e and Omega_b

homega_b <- diag(1/sqrt(diag(Sigma_b)))
hsigma_b <- diag(sqrt(diag(Sigma_b)))
res <- eigen(hsigma_b %*% Omega_e %*% hsigma_b)
D_e <- diag(res$values)
V_e <- homega_b %*% res$vectors
f2(V_e %*% D_e %*% t(V_e) - Omega_e)
f2(V_e %*% t(V_e) - Omega_b)
iV_e <- solve(V_e)
f2(iV_e %*% V_e - eye(pY))

VI <- V_e %x% eye(pX)
icov2 <- VI %*% (D_e %x% xtx + eye(pX * pY)) %*% t(VI)
f2(icov - icov2)

resX <- svd(rbind(X, zeros(pX - n, pX)))
V_x <- resX$v
D_x <- diag(resX$d)
#f2(V_x %*% D_x^2 %*% t(V_x) - xtx)
f2(V_x %*% t(V_x) - eye(pX))
VV <- V_e %x% V_x
iVV <- iV_e %x% t(V_x)
f2(iVV %*% VV - eye(pX * pY))


icov2 <- VV %*% (D_e %x% D_x^2 + eye(pX * pY)) %*% t(VV)
f2(icov2 - icov)
image(VV %*% t(VV))
post_cov2 <- t(iVV) %*% solvediag(D_e %x% D_x^2 + eye(pX * pY), iVV)
f2(post_cov2 - post_cov)


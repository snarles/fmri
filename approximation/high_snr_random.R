####
##  Randomized laplace approximation... most accurate yet??
####

source("approximation//high_snr_random_source.R")
library(magrittr)
######
### Demo
######

## small K so mc is OK
p <- 5; Sigma <- 5 * cov(randn(2*p, p))
300 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.),mc3c(Sigma,.))}

p <- 2; Sigma <- 20 * cov(randn(2*p, p))
300 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.),mc3c(Sigma,.))}

p <- 20; Sigma <- 0.5 * cov(randn(2*p, p))
100 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.))}

p <- 50; Sigma <- 0.2 * cov(randn(2*p, p))
100 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.))}


###
## larger K, don't use mc
###
p <- 10; n <- 5
Sigma <- 7 * cov(randn(2*p, p))
K <- 5000
mc2(Sigma, K)
mc3b(Sigma, K)
mc3c(Sigma, K)
## check approximations
Omega <- solve(Sigma)
Amat <- eye(p) - solve(eye(p) + Omega)
f2(eye(p) + Omega, eye(p))
f2(Amat, Omega)

###
# conditional moments of ||y-mu*||_Sigma
###
p <- 10; n <- 5; d <- p
Sigma <- 7 * cov(randn(2*p, p))
K <- 3000
p <- dim(Sigma)[1]
Omega <- solve(Sigma)
Amat <- eye(p) - solve(eye(p) + Omega)
Ha_Y <- lineId::sqrtm(eye(p) + Omega)
Ha_mu <- lineId::sqrtm(Amat)
nmS <- function(v) t(v) %*% Sigma %*% v
y <- as.numeric(Ha_Y %*% rnorm(p))
en0 <- as.numeric(t(y)%*%Amat%*%Sigma%*%Amat%*%y + TR(Amat%*%Sigma))
vn0 <- as.numeric((4*t(y)%*%Amat%*%Sigma%*%Amat%*%Sigma%*%Amat%*%y+
  2 * TR(Amat%*%Sigma%%Amat%*%Sigma)))
(ss <- (-(2 * en0/d) + sqrt((2*en0/d)^2 + 2 * vn0/d))/2)
(ll <- en0/ss - d)

n0s <- sapply(1:1e5, function(i) {
  y_mu <- as.numeric(Ha_mu %*% rnorm(p) + Amat %*% y)
  n0 <- t(y_mu) %*% Sigma %*% y_mu})
c(mean(n0s), en0)
c(var(n0s), vn0)

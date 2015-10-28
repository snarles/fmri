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
300 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3a(Sigma,.),mc3b(Sigma,.))}


p <- 2; Sigma <- 20 * cov(randn(2*p, p))
300 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.))}

p <- 20; Sigma <- 0.5 * cov(randn(2*p, p))
100 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.))}

p <- 50; Sigma <- 0.2 * cov(randn(2*p, p))
100 %>% {c(mc(Sigma,.),mc2(Sigma,.),mc3b(Sigma,.))}


###
## larger K, don't use mc
###
p <- 10; n <- 5
Sigma <- 7 * cov(randn(2*p, p))
K <- 3000
mc2(Sigma, K)
mc3a(Sigma, K)
## check approximations
Omega <- solve(Sigma)
Amat <- eye(p) - solve(eye(p) + Omega)
f2(eye(p) + Omega, eye(p))
f2(Amat, Omega)

####
## Test sampling dist code
####

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
library(glmnet, warn.conflicts = FALSE)
library(parallel)
source('utils/zattach.R')
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/bayes_reg.R')
source('eb_ident/source.R')

hyperpars <- list(n=3, pY= 4, pX=4 , W_X= 2, s_e= 1, s_b = 0.01, 
                  W_e= 2, rho_t = 0.9, L= 100, n_te= 100)
pars <- do.call(gen_params, hyperpars)
lambdas <- rep(1, hyperpars$pX)

nreps <- 100
Bs <- mclapply(1:nreps, function(i) {
  truth <- do.call(gen_data, pars)
  do.call2(samp_moments, truth, lambdas = lambdas,  matrix = TRUE, computeCov = FALSE)
}, mc.cores = 3)
truth <- do.call(gen_data, pars)


Bcov <- do.call2(samp_moments, truth, lambdas = lambdas, computeCov = TRUE)$Cov

x_star <- truth$X_te[1, ]
pres <- do.call(rbind, lapply(Bs, function(b) t(x_star) %*% b))
preCovE <- cov(pres)

pre_stuff <- do.call2(samp_predictive, truth, lambdas = lambdas)
preCov <- pre_stuff[[1]]$Cov
f2(preCovE - preCov)

Bs <- do.call(rbind, lapply(Bs, as.numeric))
BcovE <- cov(Bs)

plot(BcovE, Bcov)
f2(BcovE - Bcov)

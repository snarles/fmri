library(pracma)
library(MASS)
library(class)
library(parallel)
library(FNN)
library(lineId) ## use devtools::install('lineId')
#source('extrapolation/kay_method.R')
source('fingerprint/mi_est_pipeline.R')


K <- 30
d <- 18
p <- 100
sigma2 <- 0.3

xs <- randn(K, d)
amat <- randn(d, p)
ys <- (xs + sqrt(sigma2) * randn(K, d)) %*% amat
ys2 <- (xs + sqrt(sigma2) * randn(K, d)) %*% amat

# true MI between X and Y, I(X; Y)
rho_true <- 1/sqrt(1 + sigma2)
mi_true <- d * (-.5 * log(1 - rho_true^2))
mi_true

# scores matrix
pmat <- -pdist2(ys, ys2)
resample_misclassification(pmat, 1:K)

K <- dim(pmat)[1]
empirical_acc <- 1- resample_misclassification(pmat, 1:K)

mi_grid <- 3 * (1:300)/200
pi_ests <- piK(sqrt(2 * mi_grid), K)
mi_est <- max(mi_grid[(1-pi_ests) < empirical_acc ])
mi_est

component_mis <- sapply(1:p, function(i) mutinfo(ys[, i], ys2[, i]))
sum(component_mis)

ratio <- sum(component_mis)/mi_est
est_d <- p/ratio
est_d


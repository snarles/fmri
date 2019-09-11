library(pracma)
library(MASS)
library(class)
library(parallel)
library(lineId) ## use devtools::install('lineId')
source('extrapolation/kay_method.R')
source('fingerprint/mi_est_pipeline.R')


K <- 30
p <- 17
sigma2 <- 0.8

xs <- randn(K, p)
ys <- xs + sqrt(sigma2) * randn(K, p)
ys2 <- xs + sqrt(sigma2) * randn(K, p)

# true MI between X and Y, I(X; Y)
rho_true <- 1/sqrt(1 + sigma2)
mi_true <- p * (-.5 * log(1 - rho_true^2))
mi_true

# scores matrix
pmat <- -pdist2(ys, ys2)
resample_misclassification(pmat, 1:K)

# extrapolate until 0.5 acc
mi_est_pipeline2(pmat)

library(pracma)
library(MASS)
library(class)
library(parallel)
library(lineId) ## use devtools::install('lineId')
source('extrapolation/kay_method.R')

K <- 10
p <- 17
sigma2 <- 0.4

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
kgrid <- floor(exp((1:200)/10))
accs_kay <- kernel_extrap(pmat, kgrid)
accs_kay

kchosen <- min(kgrid[accs_kay < 0.5])
acc_at_k <- max(accs_kay[accs_kay < 0.5])
list(kchosen=kchosen, acc_at_k=acc_at_k)

## get MI estimate
mi_est <- Ihat_LI(1-acc_at_k, kchosen)
piK(mi_est/2, kchosen)
mi_est


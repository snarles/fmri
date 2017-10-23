source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")
source("par2/objective_function.R")
source("extrapolation/basis_source.R")
library(parallel)

p <- 10
sigma2 <- 0.21
sigma2_tr <- sigma2
K <- 1600
Ktarg <- 1600

repno <- 25
subfun <- function (repno) {
  set.seed(repno)
#  sigma2 <- sigma2s[repno]
#  sigma2_tr <- sigma2s[repno]
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  #pmat_sub <- -pdist2(ys[1:ksub, ], mu_hats[1:ksub, ])
  rSqs <- rowSums((ys - mu_hats)^2)
  counts <- countDistEx(mu_hats, ys, rSqs)
  accs <- count_acc(counts, Ktarg)
  accs
}

acs <- sapply(1:50, subfun)
mean(acs)
sd(acs)

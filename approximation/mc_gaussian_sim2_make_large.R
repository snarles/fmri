library(parallel)
library(abind)

mcc <- 40

source("approximation/gaussian_identity_finsam.R")
source("extrapolation/ku_source.R")

p <- 10

sigma2 <- 0.25
sigma2_tr <- sigma2
Ks <- c(16) * 1000
times <- 0 * Ks


for (ii in 1:length(times)) {
  K <- Ks[ii]
  ksub <- K/4

  t1 <- proc.time()
  
  nsplines <- 100
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- rev(1 - knts[-c(1, nsplines + 2)]^2)
  MM_2_sq <- spline1_moments(knts, 1:K)
  
  nsplines <- 100
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- knts[-c(1, nsplines + 2)]
  MM_2 <- spline1_moments(knts, 1:K)
  
  nsplines <- 1000
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- knts[-c(1, nsplines + 2)]
  MM_3 <- spline1_moments(knts, 1:K)
  
  nsplines <- 1000
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- rev(1 - knts[-c(1, nsplines + 2)]^2)
  MM_3_sq <- spline1_moments(knts, 1:K)
  
  nsplines <- 10000
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- knts[-c(1, nsplines + 2)]
  MM_4 <- spline1_moments(knts, 1:K)
  
  nsplines <- 10000
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- rev(1 - knts[-c(1, nsplines + 2)]^2)
  MM_4_sq <- spline1_moments(knts, 1:K)
  
  basis_vecs <- list(mm2 = MM_2, mm2sq =MM_2_sq, mm3 = MM_3, mm3sq = MM_3_sq, mm4 = MM_4, mm4sq = MM_4_sq)
  
  preds <- matrix(0, length(sigma2s), length(basis_vecs))
  mus <- randn(K, p)
  noise1 <- randn(K, p)
  noise2 <- randn(K, p)
  accsS <- matrix(0, length(sigma2s), K)
  ys <- mus + sqrt(sigma2) * noise1
  mu_hats <- mus + sqrt(sigma2_tr) * noise2
  pmat <- -pdist2(mu_hats, ys)
  accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
  pmat_sub <- pmat[1:ksub, 1:ksub]
  accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
  preds <- numeric(length(basis_vecs))
  for (ind in 1:length(basis_vecs)) {
    MM <- basis_vecs[[ind]]
    bt <- nnls::nnls(MM[1:ksub, ], 1- accs_sub)
    preds[ind] <- 1 - (MM[K, , drop = FALSE] %*% bt$x)
  }
  times[ii] <- (proc.time() - t1)[3]
  print(list(K = K, preds = preds, true_acc = accs[K], time = times[ii]))
}


times

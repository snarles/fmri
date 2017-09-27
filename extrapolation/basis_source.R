library(parallel)
source("extrapolation/ku_source.R")
source("extrapolation/gaussian_moments.R")



rmse_table <- function(accs_hat) {
  resids <- accs_hat - facc_rep
  resids2 <- resids^2
  biases <- t(sapply(1:50, function(ind) colMeans(resids[sigma2s==sigma2s[ind], ])))
  rmses <- t(sapply(1:50, function(ind) sqrt(colMeans(resids2[sigma2s==sigma2s[ind], ]))))
  list(biases = biases, rmses = rmses)
}

get_basis_mat <- function(max.mu, kernel_sigma, mc.reps = 1e4) {
  (n_half <- ceiling(max.mu/kernel_sigma))
  (seq_half <- seq(0, max.mu, length.out = n_half))
  (kernel_mus <- c(rev(seq_half[-1]), seq_half, Inf))
  
  Xmat <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, kref, mc.reps))
  Xtarg <- sapply(kernel_mus, function(mu) gaussian_basis(mu, kernel_sigma^2, Ktarg, mc.reps))
  list(Xmat = Xmat, Xtarg = Xtarg)
}

get_pred <- function(ind, Xmat, Xtarg) {
  nnls_fit <- nnls(Xmat, accs_subs[ind, ])
  as.numeric(Xtarg %*% nnls_fit$x)
}

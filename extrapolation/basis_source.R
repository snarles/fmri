library(parallel)
source("extrapolation/ku_source.R")
source("extrapolation/gaussian_moments.R")
library(nnls)


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

get_pred <- function(accs_sub, set1) {
  nnls_fit <- nnls(set1$Xmat, accs_sub)
  as.numeric(set1$Xtarg %*% nnls_fit$x)
}

bdwid_all_preds <- function(accs_sub, basis_sets) {
  ans <- matrix(NA, length(basis_sets), nrow(basis_sets[[1]]$Xtarg))
  ntr <- length(accs_sub)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat, accs_sub)
    ans[j, ] <- set1$Xtarg %*% nnls_fit$x
  }
  ans
}

bdwid_cv_curve <- function(accs_sub, basis_sets, cv.frac = 0.5) {
  ans <- numeric()
  ntr <- length(accs_sub)
  tr.inds <- 1:floor(ntr * cv.frac)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat[tr.inds, ], accs_sub[tr.inds])
    pred <- Xmat[-tr.inds, ] %*% nnls_fit$x
    ans[j] <- sqrt(mean((pred - accs_sub[-tr.inds])^2))
  }
  ans
}

bdwid_fit_curve <- function(accs_sub, basis_sets) {
  ans <- numeric()
  ntr <- length(accs_sub)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat, accs_sub)
    ans[j] <- sqrt(mean(nnls_fit$residuals^2))
  }
  ans
}

bdwid_cv_curve <- function(accs_sub, basis_sets, cv.frac = 0.5) {
  ans <- numeric()
  ntr <- length(accs_sub)
  tr.inds <- 1:floor(ntr * cv.frac)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat[tr.inds, ], accs_sub[tr.inds])
    pred <- Xmat[-tr.inds, ] %*% nnls_fit$x
    ans[j] <- sqrt(mean((pred - accs_sub[-tr.inds])^2))
  }
  ans
}

bdwid_fit_curve <- function(accs_sub, basis_sets) {
  ans <- numeric()
  ntr <- length(accs_sub)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat, accs_sub)
    ans[j] <- sqrt(mean(nnls_fit$residuals^2))
  }
  ans
}

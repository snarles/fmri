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

get_basis_mat <- function(max.mu, kernel_sigma, kref, Ktarg, mc.reps = 1e4) {
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

bdwid_cv_curve_f <- function(accs_sub, basis_sets, cv.frac = 0.5) {
  ans <- numeric()
  ntr <- length(accs_sub)
  tr.inds <- 1:floor(ntr * cv.frac)
  for (j in 1:length(basis_sets)) {
    set1 <- basis_sets[[j]]
    Xmat <- set1$Xmat
    nnls_fit <- nnls(Xmat[tr.inds, ], accs_sub[tr.inds])
    pred <- Xmat[nrow(Xmat), , drop = FALSE] %*% nnls_fit$x
    ans[j] <- abs(pred - accs_sub[nrow(Xmat)])
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


## basis selection code for multiple repeat

library(lineId)
cv_reg <- function(pmat, ksub2, nboot, i_chosen, basis_sets, accs = NULL) {
  ksub <- ncol(pmat)
  sub_basis_sets <- lapply(basis_sets, function(set1) {
    list(Xmat = set1$Xmat[1:ksub2, ], Xtarg = set1$Xmat[ksub, , drop = FALSE])
  })
  if (is.null(accs)) {
    accs <- 1 - resample_misclassification(pmat, i_chosen, 1:ksub)
  }
  boot_accs <- matrix(NA, nboot, ksub2)
  for (ii in 1:nboot) {
    subclasses <- sort(sample(ksub, ksub2, replace = FALSE))
    pmat_sub <- pmat[i_chosen %in% subclasses, ]
    i_chosen2 <- i_chosen[i_chosen %in% subclasses]
    remap <- numeric(ksub)
    remap[subclasses] <- 1:ksub2
    i_chosen2 <- remap[i_chosen2]
    boot_accs[ii, ] <- 1 - resample_misclassification(pmat_sub, i_chosen2, 1:ksub2)
  }
  all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
  cv_curve <- sqrt(colMeans((all_sub_preds - accs[ksub])^2))
  sel_ind <- which.min(cv_curve)
  Xmat <- basis_sets[[sel_ind]]$Xmat
  Xpred <- basis_sets[[sel_ind]]$Xtarg
  bt <- nnls::nnls(Xmat, accs)
  Xpred %*% bt$x
}

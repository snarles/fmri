#' Multivariate ridge regression
#' 
#' A \code{forward_method}.
#' Fits multivariate ridge with the same lambda for all columns
#' Computes \code{B = solve(t(X)\%*\% X + lambda * eye(p), t(X) \%*\% Y)}.
#' Uses kernel trick and parallelization.
#' @param X design matrix
#' @param Y response matrix
#' @param filt Logical vector: which responses to use
#' @param lambda L2 regularization
#' @param mc.cores number of mc cores
#' @return B coefficient matrix
#' @export
fit_ridge_kernel <- function(X, Y, filt = rep(TRUE, dim(Y)[2]),
                             lambda = 1, mc.cores = 3, ...) {
  Y <- Y[, filt]
  n <- dim(X)[1]; p <- dim(X)[2]
  gm <- paramultiply(X, t(X), mc.cores) + lambda * eye(n)
  B <- paramultiply(t(X), solve(gm, Y), mc.cores)
  B
}

#' Multivariate elastic net regression
#' 
#' A \code{forward_method}.
#' Fits elastic net ridge separately to each column
#' @import glmnet
#' @param X design matrix
#' @param Y response matrix
#' @param filt Logical vector: which responses to use
#' @param alpha Ridge (0) or Lasso (1)
#' @param rule Use \code{lambda.1se} or \code{lambda.min}
#' @param mc.cores number of mc cores
#' @return B coefficient matrix
#' @export
fit_elnet_CV <- function(X, Y, filt = rep(TRUE, dim(Y)[2]), alpha = 0,
                         rule = "lambda.1se", mc.cores = 1, ...) {
  Y <- Y[, filt]
  n <- dim(X)[1]; pX <- dim(X)[2]; pY <- dim(Y)[2]
  stuff <- function(i) {
    res <- cv.glmnet(X, Y[, i], alpha = alpha, 
                     intercept = FALSE, standardize = FALSE,
                     grouped = FALSE)
    coef.cv.glmnet(res, s = res[[rule]])[-1]
  }
  ss <- mclapply0(1:pY, stuff, mc.cores = mc.cores)
  B <- do.call(cbind, ss)
  B
}


#' Covariance estimation using residuals
#' 
#' A \code{Sigma_e_method}.
#' Estimates within-row covariances (or \code{Sigma_e}) 
#' using sample covariance of residuals,
#' applying off-diagonal shrinkage.
#' @param X design matrix
#' @param Y response matrix
#' @param B fitted coefficient matrix
#' @param filt Logical vector: which responses to use
#' @param shrink shrinkage factor
#' @param mc.cores parallelization
#' @return Estimate for covariance matrix
#' @export
residual_offdiag <- function(X, Y, B, filt = rep(TRUE, dim(Y)[2]),
                             shrink = 0.5, mc.cores = 1, ...) {
  Y <- Y[, filt]
  Yh <- paramultiply(X, B, mc.cores)
  resids <- Y - Yh
  (1- shrink) * cov(resids) + shrink * diag(diag(cov(resids)))
}

#' Produces an identity matrix
#' 
#' A \code{Sigma_t_method}.
#' @param X design matrix
#' @param Y response matrix
#' @param B fitted coefficient matrix
#' @param filt Logical vector: which responses to use
#' @param shrink shrinkage factor
#' @param mc.cores parallelization
#' @return Estimate for covariance matrix
#' @import pracma
#' @export
assume_iid <- function(X, Y, B, filt = rep(TRUE, dim(Y)[2]), ...) {
  eye(dim(X)[1])
}
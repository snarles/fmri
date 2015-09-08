#' Computing posterior likelihoods from predictive moments
#' 
#' Computes gaussian likelihoods for every (response, class) combo.
#' @param X_te Test class covariates
#' @param y_star Responses to be classified
#' @param pre_moments Predictive means and covariances for test classes
#' @param filt Logical vector: which responses to use
#' @param mc.cores Number of parallel cores
#' @return Matrix of log-likelihoods, rows = responses, columns = classes
#' @export
post_likes <- function(X_te, y_star, pre_moments,
                       filt = rep(TRUE, length(y_star)),
                       mc.cores, ...) {
  y_filt <- y_star[, filt]
  L <- dim(X_te)[1]; n_te <- dim(y_filt)[1]; pY <- dim(y_filt)[2]
  colf <- function(i) {
    Mu <- pre_moments[[i]]$Mu
    Cov <- pre_moments[[i]]$Cov
    ld <- log(det(Cov))
    resid <- t(t(y_filt) - Mu)
    ss <- solve(Cov, t(resid))
    ips <- rowSums(resid * t(ss))
    -ld - ips
  }
  res <- mclapply0(1:L, colf, mc.cores = mc.cores)
  plikes <- do.call(cbind, res)
  plikes
}



#' Predictive distribution corresponding to MLE rule
#' 
#' A \code{backward_method}.
#' MLE rule is equivalent to using the same covariance for each response
#' @param X_te Test class covariates
#' @param B fitted coefficients
#' @param Sigma_e Fitted response covariance
#' @param filt Logical vector: which responses to use
#' @return A list with moments
#' @export
pre_mle <- function(X_te, B, Sigma_e, filt, ...) {
  L <- dim(X_te)[1]
  ans <- as.list(numeric(L))
  Mus <- X_te %*% B
  for (i in 1:L) ans[[i]] <- list(Mu = Mus[i, ], Cov = Sigma_e)
  ans
}



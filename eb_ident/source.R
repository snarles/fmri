post_probs <- function(Sigma_e, Sigma_vec_B, B_hat, X, Y) {
  L <- dim(X)[1]
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  mus <- list(L)
  covs <- list(L)
  pprobs <- matrix(0, n, L)
  for (i in 1:L) {
    mus[[i]] <- X[i, , drop = FALSE] %*% B_hat
    covs[[i]] <- (diag(rep(1, p)) %x% t(X[i, ])) %*% 
      Sigma_vec_B %*% (diag(rep(1, p)) %x% t(t(X[i, ]))) + Sigma_e
    for (j in 1:n) {
      pprobs[j, i] <- -log(det(covs[[i]])) - 
        (Y[j, , drop = FALSE] - mus[[i]]) %*% solve(covs[[i]]) %*% 
        t(Y[j, , drop = FALSE] - mus[[i]])
    }  
  }
  pprobs
}
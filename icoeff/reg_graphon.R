## Empirically testing data-processing inequality for cycle densities on regular graphons

library(matpow)

## creates a regular weighted graph
gen_rwg <- function(shape = 1, rate = 1, p = 3, self.loop = TRUE, normalize = FALSE, tol = 1e-10) {
  vals <- rgamma(p^2, shape, rate)
  mat <- matrix(vals, p, p)
  if (!self.loop) diag(mat) <- 0
  while (max(abs(rowSums(mat) -1)) > tol || max(abs(colSums(mat) - 1)) > tol) {
    mat <- mat/rowSums(mat)
    mat <- t(t(mat)/colSums(mat))
  }
  if (normalize) mat <- mat/p
  # list(rowSums(mat), colSums(mat))
  mat
}

tr_power <- function(mat, pow = 3) {
  sum(diag(matpow(mat, pow, squaring = TRUE)$prod1))
}

pspec <- function(mat, max.pow = 10) {
  evs <- eigen(mat)$values
  lala <- t(t(rep(1,max.pow))) %*% evs
  lala <- lala ^ row(lala)
  Re(rowSums(lala))
}

check_reg <- function(mat, tol = 1e-10) {
  stopifnot(max(abs(rowSums(mat) -  1)) < tol)
  stopifnot(max(abs(colSums(mat) -  1)) < tol)
}

array_colsum <- function(arr) {
  ad <- dim(arr)
  mat <- matrix(arr, ad[1], ad[2] * ad[3])
  matrix(colSums(mat), ad[2], ad[3])
}

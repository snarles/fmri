library(Matrix)
library(pracma)

## finds a matrix X such that X'=AX+B
transpose_solve = function(A, B) {
  p <- dim(A)[1]
  bv <- -as(t(t(as.numeric(B))), "dgCMatrix")
  AM <- as(eye(p) %x% A, "dgCMatrix")
  inds <- (rep(1:p, p)-1)*p + (rep(1:p, each = p) - 1) + 1
  TM <- sparseMatrix(1:p^2, inds, x=rep(1, p^2))
  M <- AM - TM
  xv <- solve(M, bv)
  X <- matrix(xv, p, p)
  stopifnot
    (
      max(abs(t(X) - A %*% X - B)) < 1e-8
    )
  X
}

# B <- randn(5)
# A <- randn(5)
# X <- transpose_solve(A, B)
# t(X)
# A %*% X + B

library(Matrix)
library(pracma)

TR <- function(x) sum(diag(x))

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

## minimizes O
##  tr(OO'H)-2logdet(O)-2tr(CO)
jmle_sub_1 <- function(H, C) {
  objf <- function(O) TR(O%*%t(O)%*%H) - 2*log(abs(det(O))) - 2*TR(C %*% O)
  O1 <- transpose_solve(H, -t(C))
  O2 <- transpose_solve(-H, t(C))
  scores <- c(objf(O1), objf(O2))
  list(O1, O1)[[order(scores)[1]]]
}

## minimizes G:
##  tr(GG'Z) - tr(GC) + lambda * tr((I-G'G)^2)
jmle_sub_2a <- function(Z, C, lambda)
  transpose_solve(-(Z/lambda + eye(dim(Z)[1])), C/lambda)

## minimizes G subject to G'G = I
##  tr(GG'Z) - tr(GC)
jmle_sub_2 <- function(Z, C) {
  lambda <- 1
  G <- jmle_sub_2a(Z, C, lambda)
}

## minimizes G, A, B subject to G'G = I:
##  tr((Y - XGA)(A'A-)(Y - XGA)')
##   + tr((W - XGB)(B'B-)(W - XGB)')
##   + logdet(A'A) + logdet(B'B)
jmle <- function(X, Y, W, n.its = 20) {
  q <- dim(Y)[2]; p <- dim(X)[2]
  G <- svd(randn(p, q))$u;
  A <- randn(q); B <- randn(q)
  objf <- function(G, A, B) {
    Yh <- X %*% G %*% A; Wh <- X %*% G %*% B
    s1 <- (Y-Yh) %*% solve(t(A) %*% A, t(Y-Yh))
    s2 <- (W-Wh) %*% solve(t(B) %*% B, t(W-Wh))
    TR(s1) + TR(s2) + log(det(t(A) %*% A)) + 2 * log(det(t(B) %*% B))
  }
  for (it in 1:n.its) {
    objf(G, A, B)
    ## Optimize A
    H <- t(Y) %*% Y; C <- t(t(Y) %*% X %*% G)
    A <- solve(jmle_sub_1(H, C))
    objf(G, A, B)
    ## Optimize B
    H <- t(W) %*% W; C <- t(t(W) %*% X %*% G)
    B <- solve(jmle_sub_1(H, C))
    objf(G, A, B)
    ## Optimize G
    Z <- t(X) %*% X %*% G
    C <- (t(X) %*% t(solve(t(A), t(Y)) + solve(t(B), t(Y))))/2
    
  }
}  

n <- 50; p <- 30; q <- 20
X <- randn(n, p)
G0 <- svd(randn(p, q))$u
A0 <- randn(q, q); B0 <- randn(q, q)
Y <- (X %*% G0  + randn(n, q)) %*% A0
W <- (X %*% G0  + randn(n, q)) %*% B0

library(Matrix)
library(pracma)
library(lineId)
library(magrittr)

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
##  tr(OO'H)-2nlogdet(O)-2tr(CO)
jmle_sub_1 <- function(H, C, n) {
  objf <- function(O) TR(O%*%t(O)%*%H) - 2*n*log(abs(det(O))) - 2*TR(C %*% O)
  O1 <- transpose_solve(H/n, -t(C)/n)
  O2 <- transpose_solve(-H/n, t(C)/n)
  scores <- c(objf(O1), objf(O2))
  O <- list(O1, O2)[[order(scores)[1]]]
}

objf <- function(O) TR(O%*%t(O)%*%H) - 2*n*log(abs(det(O))) - 2*TR(C %*% O)
O1 <- O
ss <- 1e-6
objf(O1)

for (i in 1:100) {
  GOa <- 2 * H %*% O1 - 2 *n * t(O1) -2*t(C)
  GOb <- 2 * H %*% O1 + 2 *n * t(O1) -2*t(C)
  s1 <- objf(O1 - ss * GOa)
  s2 <- objf(O1 - ss * GOb)
  O1 <- O1 - ss * list(GOa,GOb)[[order(c(s1,s2))[1]]]
  print(objf(O1))  
}

## minimizes G, A, B subject to G'G = I:
##  tr((Y - XGA)(A'A-)(Y - XGA)')
##   + tr((W - XGB)(B'B-)(W - XGB)')
##   + logdet(A'A) + logdet(B'B)
jmle <- function(X, Y, W, n.its = 20) {
  q <- dim(Y)[2]; p <- dim(X)[2]; n <- dim(X)[1]
  G <- randn(p, q);
  A <- randn(q); B <- randn(q)
  objf <- function(G, A, B) {
    Yh <- X %*% G %*% A; Wh <- X %*% G %*% B
    s1 <- (Y-Yh) %*% solve(t(A) %*% A, t(Y-Yh))
    s2 <- (W-Wh) %*% solve(t(B) %*% B, t(W-Wh))
    TR(s1) + TR(s2) + n*log(det(t(A) %*% A)) + n*log(det(t(B) %*% B))
#     TR(solve(t(A) %*% A, t(Y) %*% Y)) +
#       TR(solve(t(B) %*% B, t(W) %*% W)) -
#       2*TR(G %*% solve(t(A), t(Y)) %*% X + G%*% solve(t(B), t(W)) %*%X) +
#       2*TR(G %*% t(G) %*% t(X) %*% X) +
#       n* log(det(t(A) %*% A)) + n * log(det(t(B) %*% B))
  }
  for (it in 1:n.its) {
    objf(G, A, B)
    ## Optimize A
    H <- t(Y) %*% Y; C <- t(t(Y) %*% X %*% G); O <- solve(A)
    A <- solve(jmle_sub_1(H, C, n))
    objf(G, A, B)
    ## Optimize B
    H <- t(W) %*% W; C <- t(t(W) %*% X %*% G)
    B <- solve(jmle_sub_1(H, C, n))
    objf(G, A, B)
    ## Optimize G
    Z <- t(X) %*% X
    C <- t(X) %*% t(solve(t(A), t(Y)) + solve(t(B), t(W)))/2
    G <- solve(Z, C)
    objf(G, A, B)
  }
  list(G = G, A = A, B= B)
}


## TEST

n <- 1000; p <- 2; q <- 2
X <- randn(n, p)
G0 <- randn(p, q)
A0 <- randn(q, q); B0 <- randn(q, q)
sigma <- 1
# MA0 <- G0 %*% A0 %*% solve(t(A0) %*% A0, t(G0 %*% A0))
# MB0 <- G0 %*% B0 %*% solve(t(B0) %*% B0, t(G0 %*% B0))
# f2(MA0, MB0)
# f2(MA0, G0 %*% t(G0))

Y <- (X %*% G0  + sigma * randn(n, q)) %*% A0
W <- (X %*% G0  + sigma * randn(n, q)) %*% B0
sol <- jmle(X, Y, W, 20)
f2(G0 %*% t(G0), sol$G %*% t(sol$G)) # not identifiable!
#f2(G0 %*% t(G0))
f2(G0 %*% A0, sol$G %*% sol$A)
#f2(G0 %*% A0)
f2(G0 %*% B0, sol$G %*% sol$B)
#f2(G0 %*% B0)

f2(A0)
f2(sol$A)


objf <- function(G, A, B) {
  Yh <- X %*% G %*% A; Wh <- X %*% G %*% B
  s1 <- (Y-Yh) %*% solve(t(A) %*% A, t(Y-Yh))
  s2 <- (W-Wh) %*% solve(t(B) %*% B, t(W-Wh))
  TR(s1) + TR(s2) + n*log(det(t(A) %*% A)) + n*log(det(t(B) %*% B))
}

truth <- list(G=G0, A=A0, B=B0)
objf(G0, A0, B0)
sol %$% objf(G, A, B)
sol %$% log(det(t(A) %*% A))
truth %$% log(det(t(A) %*% A))

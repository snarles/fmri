library(pracma)
library(lineId)
library(magrittr)
source("topological//rotmax.R")

TR <- function(x) sum(diag(x))

init_est <- function(X, Y, W) {
  q <- dim(Y)[2]; p <- dim(X)[2]; n <- dim(X)[1]
  #A <- randn(q); B <- randn(q);
  PX <- X %*% solve(t(X) %*% X, t(X))
  cX <- solve(t(X) %*% X, t(X))
  B_Y <- cX %*% Y; Yh <- PX %*% Y
  B_W <- cX %*% W; Wh <- PX %*% W
  SigmaY <- cov(Y - Yh)
  SigmaW <- cov(W - Wh)
  A <- sqrtm(SigmaY)
  gg <- rotmin(B_Y%*%isqrtm(SigmaY), B_W%*%isqrtm(SigmaW))
  B <- gg %*% sqrtm(SigmaW)
  M_Y <- B_Y %*% solve(SigmaY, t(B_Y))
  M_W <- B_W %*% solve(SigmaW, t(B_W))
  list(A=A, B=B, M_Y = M_Y, M_W = M_W) 
}

sep_mle_of <- function(X, Y, W) {
  q <- dim(Y)[2]; p <- dim(X)[2]; n <- dim(X)[1]
  #A <- randn(q); B <- randn(q);
  PX <- X %*% solve(t(X) %*% X, t(X))
  cX <- solve(t(X) %*% X, t(X))
  B_Y <- cX %*% Y; Yh <- PX %*% Y
  B_W <- cX %*% W; Wh <- PX %*% W
  SigmaY <- cov(Y - Yh)
  SigmaW <- cov(W - Wh)
  TR((Y-Yh) %*% solve(SigmaY, t(Y-Yh))) + 
    TR((W-Wh) %*% solve(SigmaW, t(W-Wh))) +
    n*log(det(SigmaY)) + n*log(det(SigmaW))
}

jmle <- function(X, Y, W, n.its = 200, ss = 1, ssm = 0.5,
                 init = init_est(X, Y, W)) {
  q <- dim(Y)[2]; p <- dim(X)[2]; n <- dim(X)[1]
  #A <- randn(q); B <- randn(q);
  PX <- X %*% solve(t(X) %*% X, t(X))
  A <- init$A; B <- init$B
  Ai <- solve(A); Bi <- solve(B)
  YtY <- t(Y) %*% Y
  WtW <- t(W) %*% W
  YxY <- t(Y) %*% PX %*% Y
  WxW <- t(W) %*% PX %*% W
  YxW <- t(Y) %*% PX %*% W
  objf <- function(A, B) {
    TR(solve(t(A) %*% A, t(Y) %*%Y) + solve(t(B) %*%B, t(W) %*% W)) -
            .5 * TR(PX %*% (Y %*% solve(A) + W %*% solve(B)) %*% 
                      t(Y %*% solve(A) + W %*% solve(B))) + 
            n*log(det(t(A) %*% A)) + n*log(det(t(B) %*% B))
  }
  objf_i <- function(Ai, Bi) objf(solve(Ai), solve(Bi))
  gradAi <- function(Ai, Bi)
    2 * YtY %*% Ai - YxY %*% Ai - YxW %*% Bi - 2 * n*t(solve(Ai))
  gradBi <- function(Ai, Bi)
    2 * WtW %*% Bi - WxW %*% Bi - t(YxW) %*% Ai - 2 * n*t(solve(Bi))

  for (i in 1:n.its) {
    old_obj <- objf_i(Ai, Bi)
    gAi <- gradAi(Ai, Bi)
    gBi <- gradBi(Ai, Bi)
    new_obj <- objf_i(Ai-ss*gAi, Bi-ss*gBi)
    if (new_obj > old_obj) {
      ss <- ss * ssm
    } else {
      Ai <- Ai-ss*gAi; Bi <- Bi-ss*gBi      
      #if (i %% 10==0) print(new_obj)
    }
  }
  
  A <- solve(Ai); B <- solve(Bi)
  G <- 1/2 * solve(t(X) %*%X, t(X)) %*% 
    (Y %*% solve(A) + W %*% solve(B))
  list(G = G, A = A, B= B, objf = objf)
}


objf0 <- function(G, A, B) {
  Yh <- X %*% G %*% A; Wh <- X %*% G %*% B
  s1 <- (Y-Yh) %*% solve(t(A) %*% A, t(Y-Yh))
  s2 <- (W-Wh) %*% solve(t(B) %*% B, t(W-Wh))
  TR(s1) + TR(s2) + n*log(det(t(A) %*% A)) + n*log(det(t(B) %*% B))
}
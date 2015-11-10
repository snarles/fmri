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
  list(A=A, B=B) 
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
      if (i %% 10==0) print(new_obj)
    }
  }
  
  A <- solve(Ai); B <- solve(Bi)
  G <- 1/2 * solve(t(X) %*%X, t(X)) %*% 
    (Y %*% solve(A) + W %*% solve(B))
  list(G = G, A = A, B= B, objf = objf)
}



## TEST

n <- 1000; p <- 30; q <- 20
X <- randn(n, p)
G0 <- randn(p, q)
A0 <- randn(q, q); B0 <- randn(q, q)
sigma <- 1
truth <- list(G=G0, A=A0, B=B0)

# MA0 <- G0 %*% A0 %*% solve(t(A0) %*% A0, t(G0 %*% A0))
# MB0 <- G0 %*% B0 %*% solve(t(B0) %*% B0, t(G0 %*% B0))
# f2(MA0, MB0)
# f2(MA0, G0 %*% t(G0))

Y <- (X %*% G0  + sigma * randn(n, q)) %*% A0
W <- (X %*% G0  + sigma * randn(n, q)) %*% B0
sol <- jmle(X, Y, W, 100)
sol0 <- jmle(X, Y, W, 100, init=truth)
sep <- sep_mle_of(X, Y, W)

#f2(G0 %*% t(G0), sol0$G %*% t(sol0$G))
f2(G0 %*% t(G0), sol$G %*% t(sol$G))
#f2(G0 %*% t(G0))
f2(G0 %*% A0, sol$G %*% sol$A)
#f2(G0 %*% A0)
f2(G0 %*% B0, sol$G %*% sol$B)
#f2(G0 %*% B0)

f2(A0)
f2(sol$A)


objf0 <- function(G, A, B) {
  Yh <- X %*% G %*% A; Wh <- X %*% G %*% B
  s1 <- (Y-Yh) %*% solve(t(A) %*% A, t(Y-Yh))
  s2 <- (W-Wh) %*% solve(t(B) %*% B, t(W-Wh))
  TR(s1) + TR(s2) + n*log(det(t(A) %*% A)) + n*log(det(t(B) %*% B))
}

truth %$% objf0(G, A, B)
(jof <- sol %$% objf0(G, A, B))
sol0 %$% objf0(G, A, B)

(lr <- 1/2*(jof - sep))

sol %$% log(det(t(A) %*% A))
truth %$% log(det(t(A) %*% A))

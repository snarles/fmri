source("topological/large_n.R")


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

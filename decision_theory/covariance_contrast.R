####
##  Estimating Sigma in a mulvariate regression
####
library(lineId);library(pracma)
library(MASS);library(magrittr)

TR <- function(a) sum(diag(a))
xi_mat <- function(X, varw) {
  G <- X %*% t(X)
  res <- eigen(G)
  ws <- res$values; V <- res$vectors
  wwt <- function(gm) (sum(pmax(0, (gm-ws)/varw))-1)^2
  or <- optimise(wwt, c(0, varw * max(ws)))
  gm <- or$minimum
  V %*% diag(pmax(0, (gm-ws)/varw)) %*% t(V)
}

## offdiagonal shrinkage
ofs <- function(A, shrink) {
  (1-shrink) * A + shrink * diag(diag(A))
}

cor_dist <- function(A, B) {
  f2(A/sqrt(f2(A)) - B/sqrt(f2(B)))
}

cor_dist_ofs <- function(A, Ahat, shrinks) {
  sapply(shrinks, function(s) cor_dist(A, ofs(Ahat, s)))
}

n2A <- function(x, A, num = FALSE) {
  ans <- t(x) %*% A %*% x
  if (num) ans <- as.numeric(num)
  ans
}

####
##   Random X, no repeats
####

params0 <- list(n = 30, pY=60, pX = 60, W_X=3, s_e=3,
                s_b=300, df_b=10, W_e=2)
zattach(params0)
params <- do.call(gen_params, params0)
dat <- do.call(gen_data, params)
zattach(dat)
E <- Y - X %*% B
Bhat <- fit_ridge_kernel(X, Y, lambda=50)
f2(B - Bhat)
Ehat <- Y - X %*% Bhat
SigmaY <- t(Y) %*% Y/n
Sigmahat <- residual_offdiag(X, Y, Bhat, shrink=0)
Sigma_orc <- t(E) %*% E/n


# f2(Sigma_e)
# f2(SigmaY, Sigma_e)
# f2(Sigmahat, Sigma_e)
# f2(Sigmahat2, Sigma_e)
# f2(Sigma_orc, Sigma_e)



cor_dist(SigmaY, Sigma_e)
cor_dist(ofs(SigmaY, .5), Sigma_e)
cor_dist(ofs(SigmaY, 1), Sigma_e)

cor_dist(Sigmahat, Sigma_e)
cor_dist(ofs(Sigmahat, .5), Sigma_e)
cor_dist(ofs(Sigmahat, 1), Sigma_e)

xi <- xi_mat(X, 50000)
Sigmahat2 <- t(Ehat) %*% xi %*% Ehat
sum(xi * (X %*% t(X)))
TR(xi)
TR(xi%*%xi)

Sigma_xi_E <- t(E) %*% xi %*% E

cor_dist(Sigmahat2, Sigma_e)
cor_dist(ofs(Sigmahat2, .5), Sigma_e)
cor_dist(ofs(Sigmahat2, 1), Sigma_e)

cor_dist(Sigma_orc, Sigma_e)
cor_dist(ofs(Sigma_orc,.5), Sigma_e)
cor_dist(ofs(Sigma_orc, 1), Sigma_e)

cor_dist(Sigma_xi_E, Sigma_e)
cor_dist(ofs(Sigma_xi_E, .5), Sigma_e)
cor_dist(ofs(Sigma_xi_E, 1), Sigma_e)

####
##  X with near-repeats
####

n_h <- 20; pX <- 100; pY <- 50; n <- 2 * n_h
X_h <- randn(n_h, pX)
X <- rbind(X_h, X_h) + 0.1 * randn(n, pX)
B <- randn(pX, pY)
Sigma_e <- cov(randn(2 * pY, pY))
E <- mvrnorm(n, rep(0, pY), Sigma_e)
Y <- X %*% B + E
avgm <- eye(n)/n
diffm <- matrix(c(1,-1,-1,1), 2, 2) %x% eye(n_h)/n
Xi <- xi_mat(X, 100)
image(Xi)
Sa <- n2A(Y, avgm)
Sb <- n2A(Y, diffm)
Sc <- n2A(Y, Xi)

cor_dist_ofs(Sigma_e, Sa, 0:10/10) %>% round(digits = 3)
cor_dist_ofs(Sigma_e, Sb, 0:10/10) %>% round(digits = 3)
cor_dist_ofs(Sigma_e, Sc, 0:10/10) %>% round(digits = 3)

####
##  X with clusters
####

pX <- 100; pY <- 50; n <- 40; k <- 10
X0 <- randn(k, pX)
X <- X0[sample(k, n, TRUE), ] + 0.1 * randn(n, pX)
B <- randn(pX, pY)
Sigma_e <- cov(randn(2 * pY, pY))
E <- mvrnorm(n, rep(0, pY), Sigma_e)
Y <- X %*% B + E
avgm <- eye(n)/n
Xi <- xi_mat(X, 100)
#image(Xi)
Bhat <- lineId::fit_ridge_kernel(X, Y, lambda = 0.1)
resid <- Y - X %*% Bhat
SaY <- n2A(Y, avgm)
SaR <- n2A(resid, avgm)
SbY <- n2A(Y, Xi)
SbR <- n2A(resid, Xi)

cor_dist_ofs(Sigma_e, SaY, 0:10/10) %>% round(digits = 3)
cor_dist_ofs(Sigma_e, SaR, 0:10/10) %>% round(digits = 3)
cor_dist_ofs(Sigma_e, SbY, 0:10/10) %>% round(digits = 3)
cor_dist_ofs(Sigma_e, SbR, 0:10/10) %>% round(digits = 3)


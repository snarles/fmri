####
##  Estimating Sigma in a mulvariate regression
####

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

library(lineId)
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

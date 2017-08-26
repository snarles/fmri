library(pracma)
library(MASS)
library(magrittr)
library(parallel)
f2 <- function(x, y = 0) sum((x- y)^2)
zattach <- function(ll) {
  for (i in 1:length(ll)) {
    assign(names(ll)[i], ll[[i]], envir=globalenv())
  }
}

compute_st <- function(vs, gm, ws, ts) {
  lv <- length(vs)
  lt <- length(ts)
  vv <- repmat(vs, lt, 1)
  tt <- repmat(t(t(ts)), 1, lv)
  wt <- repmat(t(t(ws * ts)), 1, lv)
  ww <- repmat(t(t(ws)), 1, lv)
  temp <- wt/(1 + tt * vv)
  zs <- -1/vs + gm * colSums(temp)
  ms <- vs/gm + ((1/gm) - 1)/zs
  vp <-   1/(1/(vs^2) - gm * colSums(ww * tt^2/(1 + tt * vv)^2))
  list(vs = vs, zs = zs, vp = vp, ms = ms)
}


simulate_ridge_risk_cov <- function(Sigma, n, p, alpha, lambda) {
  X <- mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
  w <- alpha * rnorm(p)/sqrt(p)
  y <- X %*% w + rnorm(n)
  ip <- t(X) %*% y
  what <- solve(t(X) %*% X + n * lambda * eye(p), ip)
  Xte <- mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
  yte <- Xte %*% w + rnorm(n)
  yh <- Xte %*% what
  ytr <- X %*% w + rnorm(n)
  yh_tr <- X %*% what
  c(f2(yte, yh)/n, f2(ytr, yh_tr)/n)
}

risk_formula <- function(ts, ws, alpha2, gamma, n = 1e3, vs = linspace(0, 0.5, 10000)) {
  alpha <- sqrt(alpha2); gm <- gamma
  grid_size <- 1e3
  p <- n * gm
  res <- compute_st(vs, gm, ws, ts); zs <- res$zs; vp <- res$vp; ms <- res$ms
  ls <- -zs
  pred_risk <- (1 + (ls * alpha^2/gm - 1) * (1 - ls * vp/vs))/(ls * vs)
  list(lambdas = ls, risk = pred_risk)
}


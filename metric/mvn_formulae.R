library(pracma); library(MASS); library(lineId)

TR <- function(A) sum(diag(A))

## moments of X'Y for (X, Y) ~ N(mu, Sigma)
inner_product_moments <- function(Sigma, mu = rep(0, dim(Sigma)[1]), naive = FALSE) {
  p.5 <- dim(Sigma)[1]/2
  if (naive) {
    samp <- mvrnorm(1e6, mu=mu, Sigma = Sigma)
    prods <- rowSums(samp[, 1:p.5] * samp[, -(1:p.5)])
    return(
      list(m = mean(prods), m2 = mean(prods^2), v = var(prods))
      )
  }
  muX <- mu[1:p.5]; muY <- mu[-(1:p.5)]
  S_X <- Sigma[1:p.5, 1:p.5]
  S_XY <- Sigma[1:p.5, -(1:p.5)]; S_YX <- t(S_XY)
  S_Y <- Sigma[-(1:p.5), -(1:p.5)]
  Si_X <- solve(S_X)
  S_Y_X <- S_Y - S_YX %*% solve(S_X, S_XY)
  m <- TR(S_XY) + sum(muX * muY)
  m2 <- TR(S_XY %*% S_XY) + TR(S_XY)^2 + TR(S_X %*% S_Y) +
    sum(muX * muY)^2 + 2 * sum(muX * muY) * TR(S_XY) +
    t(muX) %*% S_Y %*% muX + t(muY) %*% S_X %*% muY +
    2 * t(muX) %*% S_YX %*% muY
  v <- TR(S_XY %*% S_XY) + TR(S_X %*% S_Y) +
    t(muX) %*% S_Y %*% muX + t(muY) %*% S_X %*% muY +
    2 * t(muX) %*% S_YX %*% muY
  list(m = m[1], m2 = m2[1], v = v[1])
}

## E[X'Y X'AX], E[X'Y X'AY] centered
more_inner_product_moments <- function(Sigma, A = eye(dim(Sigma)[1]), naive = FALSE) {
  p.5 <- dim(Sigma)[1]/2
  if (naive) {
    samp <- mvrnorm(1e6, mu=rep(0, 2 * p.5), Sigma = Sigma)
    Xs <- samp[, 1:p.5]
    Ys <- samp[, -(1:p.5)]
    xty_s <- rowSums(Xs * Ys)
    xtAy_s <- rowSums(Xs * (Ys %*% t(A)))
    xtAx_s <- rowSums(Xs * (Xs %*% t(A)))
    return(
      list(xty_xtAy = mean(xty_s * xtAy_s),
           xty_xtAx = mean(xty_s * xtAx_s))
    )
  }
  S_X <- Sigma[1:p.5, 1:p.5]
  S_XY <- Sigma[1:p.5, -(1:p.5)]; S_YX <- t(S_XY)
  S_Y <- Sigma[-(1:p.5), -(1:p.5)]
  Si_X <- solve(S_X)
  S_Y_X <- S_Y - S_YX %*% solve(S_X, S_XY)
  B_YX <- S_YX %*% solve(S_X)
  xty_xtAx <- TR(S_YX %*% (A + t(A)) %*% S_X) +
    TR(S_XY) * TR(A %*% S_X)
  xty_xtAy <- TR(A %*% S_YX %*% (B_YX + t(B_YX)) %*% S_X) +
    TR(S_YX) * TR(A %*% S_YX) +
    TR(S_X %*% S_Y_X %*% t(A))
  list(xty_xtAy = xty_xtAy, xty_xtAx = xty_xtAx)
}

####
##  TESTS
####

p <- 2 * 5
mu <- rep(0, p)
mu <- rnorm(p)
Sigma <- cov(randn(2*p, p))
A <- randn(p/2)
inner_product_moments(Sigma, mu, naive = TRUE)
inner_product_moments(Sigma, mu)
more_inner_product_moments(Sigma, A, naive = TRUE)
more_inner_product_moments(Sigma, A)


library(pracma); library(MASS); library(lineId)

TR <- function(A) sum(diag(A))
TR2 <- function(a) sum(diag(a %*% a))

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
  xty_xtAy <- TR(S_YX) * TR(A %*% S_YX) +
    TR(A %*% S_YX %*% S_YX) +
    TR(A %*% S_Y %*% S_X)
  list(xty_xtAy = xty_xtAy, xty_xtAx = xty_xtAx)
}


## moments of ||x||^2, ||y||^2, and (x, y) ~ N(0, Sigma)
norms_inner_products <- function(Sigma, naive = FALSE) {
  p.5 <- dim(Sigma)[1]/2
  mu <- rep(0, p.5 * 2)
  if (naive) {
    samp <- mvrnorm(1e6, mu=mu, Sigma = Sigma)
    xs <- samp[, 1:p.5]
    ys <- samp[, -(1:p.5)]
    nmsx <- rowSums(xs^2)
    nmsy <- rowSums(ys^2)
    nms <- cbind(nmsx, nmsy)
    return(
      list(m = colMeans(nms), v = cov(nms))
    )
  }
  S_X <- Sigma[1:p.5, 1:p.5]
  S_XY <- Sigma[1:p.5, -(1:p.5)]; S_YX <- t(S_XY)
  S_Y <- Sigma[-(1:p.5), -(1:p.5)]
  Si_X <- solve(S_X)
  S_Y_X <- S_Y - S_YX %*% solve(S_X, S_XY)
  m <- c(TR(S_X), TR(S_Y))
  v <- 2* cbind(c(TR2(S_X), TR(S_XY %*% S_YX)),
             c(TR(S_XY %*% S_YX), TR2(S_Y)))
  list(m = m, v = v)
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
norms_inner_products(Sigma, naive = TRUE)
norms_inner_products(Sigma)

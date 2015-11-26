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

####
##  TESTS
####

p <- 2 * 10
mu <- rep(0, p)
mu <- rnorm(p)
Sigma <- cov(randn(2*p, p))
inner_product_moments(Sigma, mu, naive = TRUE)
inner_product_moments(Sigma, mu)


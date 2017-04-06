## compute gaussian overlaps

library(pracma)
library(AlgDesign)
library(MASS)

mvdnorm <- function(x, mu, sigma2) {
  p <- length(mu)
  d2 <- rowSums((x - mu)^2)
  (1/sqrt(2 * pi * sigma2))^p * exp(-d2/(2 * sigma2))
}

max_d <- function(x, mus, sigma2) {
  mus2 <- split(t(mus), rep(1:nrow(mus), each = ncol(mus)))
  ds <- lapply(mus2, function(mu) mvdnorm(x, mu, sigma2))
  do.call(pmax, ds)
}

X <- as.matrix(gen.factorial(c(2,2,2)))

B1 <- matrix(0, 3, 3); B1[1,] <- 1
B2 <- eye(3)

Y1 <- X %*% B1
Y2 <- X %*% B2

res <- 91
pts <- as.matrix(gen.factorial(rep(res, 3)))/res * 7
(delta <- pts[2,1] - pts[1,1])

lala <- mvdnorm(pts, rep(0, 3), 1)
delta^3 * sum(lala)

sigmas <- 1:100/10

ans <- t(sapply(sigmas, function(s) {
  c(s, 
   delta^3 * sum(max_d(pts, Y1, s))/8,
   delta^3 * sum(max_d(pts, Y2, s))/8)
}))

matplot(ans[, 1], ans[, 2:3], ylim = c(0, 1))

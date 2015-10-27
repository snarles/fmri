library(pracma)

Kfunc <- function(Sigma) {
  d <- dim(Sigma)[1]
  1/sqrt(det(Sigma)) * 2*pi^(d/2)/d/gamma(d/2)
}

true_q <- function(Sigma, y, pp, L, mc.its = 1e5) {
  lprob <- -log(1-pp)/L
  yy <- randn(p, mc.its) + y
  nms <- colSums((pracma::sqrtm(Sigma)$B %*% yy)^2)
  sort(nms)[floor(lprob * mc.its)]
}

guess_q <- function(Sigma, y, pp, L) {
  Ks <- Kfunc(Sigma)
  lprob <- -log(1-pp)/L
  yn <- sqrt(sum(y^2))
  (lprob/((2 * pi)^(-d/2) * exp(-yn^2/2))/Ks)^(2/d)
}

## finds upper and lower bounds for x_pp
find_x_bar <- function(Sigma, y, pp, L, reso = 1e5, fudge = 10) {
  Ks <- Kfunc(Sigma)
  lprob <- -log(1-pp)/L
  p <- dim(Sigma)[1]
  res <- eigen(Sigma)
  xi <- t(res$vectors) %*% y
  gs <- res$values
  yn <- sqrt(sum(y^2))
  ## UPPER BOUND
  ## decide search range for ls
  gqs <- guess_q(Sigma, y, pp, L)
  lprox <- sqrt(sum(xi^2/gs)/min(gqs)/fudge) + 1/min(gs)
  ls <- seq(0, 2 * lprox, length.out = reso)
  lmat <- repmat(ls, p, 1)
  gmat <- repmat(t(t(gs)), 1, length(ls))
  ximat <- repmat(xi, 1, length(ls))
  xs <- colSums(ximat^2 * gmat/(lmat * gmat + 1)^2)
  stopifnot(sum(is.na(xs))==0)
  yns <- sqrt(colSums(ximat^2 * (1 - (1/(1 + lmat * gmat)))^2))
  pbs <- Ks * xs^(d/2) * (2 * pi)^(-d/2) * exp(-yns^2/2)
  #plot(ls, xs)
  #plot(xs, yns)
  #plot(xs, log(pbs))
  xs_aug <- c(xs, rep(NA, length(lprob)))
  pbs_aug <- c(pbs, lprob)
  xso <- xs_aug[order(pbs_aug)]
  pbso <- pbs_aug[order(pbs_aug)]
  na.inds <- which(is.na(xso))
  xso[na.inds] <- (xso[na.inds - 1] + xso[na.inds + 1])/2
  stopifnot(sum(is.na(xso)) == 0)
  temp <- xso[na.inds]
  temp2 <- pbso[na.inds]
  xdown <- temp[match(lprob, temp2)]
  list(xdown = xdown)
}


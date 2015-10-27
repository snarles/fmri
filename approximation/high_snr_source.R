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
find_x_bar <- function(Sigma, y, pp, L, res = reso, fudge = 10) {
  Ks <- Kfunc(Sigma)
  lprob <- -log(1-pp)/L
  p <- dim(Sigma)[1]
  res <- eigen(Sigma)
  xi <- t(res$vectors) %*% y
  gs <- res$values
  yn <- sqrt(sum(y^2))
  ## initial guess
  gqs <- guess_q(Sigma, y, pp, L)
  lprox <- sqrt(sum(xi^2/gs)/min(gqs)/fudge) + 1/min(gs)
  ls <- seq(-2 * lprox, 2 * lprox, length.out = reso)
  lmat <- repmat(ls, p, 1)
  gmat <- repmat(t(t(gs)), 1, length(ls))
  ximat <- repmat(xi, 1, length(ls))
  xs <- colSums(ximat^2 * gmat/(lmat * gmat + 1)^2)
  yns <- sqrt(colSums(ximat^2 * (1 - (1/(1 + lmat * gmat)))^2))
  10 %>% {plot(xs[xs < .], yns[xs < .])}
  #plot(ls, xs, ylim = c(0, 2*max(gqs)), type = "l")
  diff <- xs[-1]-xs[-length(xs)]
  cut1 <- ls[min(which(diff < 0))]
  cut2 <- ls[max(which(diff > 0))]
  xs1 <- xs[ls < cut1]
  yns1 <- yns[ls < cut1]
  xs2 <- xs[ls > cut2]
  yns2 <- yns[ls > cut2]
  plot(xs1, yns1, type = "o")
  plot(xs2, yns2, type = "o")
  
  plot(ls, xs * ((1:length(xs) < cut1) + (1:length(xs) > cut2)), type = "l")
  
}


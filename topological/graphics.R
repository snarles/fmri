library(lineId)
library(magrittr)

## grf on points with coordinates x, cov = exp(-lambda d(x_i, x_j)^2/2)
grf_grid <- function(xs, lambda = 2, k.basis = 5) {
  kseq <- (-k.basis):k.basis
  kgrid <- cbind(rep(kseq, each=(2*k.basis + 1)), rep(kseq, (2*k.basis + 1)))
  klen <- dim(kgrid)[1]
  zs <- rnorm(klen) + 1i * rnorm(klen)
  zs <- zs * exp(-lambda/4 * rowSums(kgrid^2))
  mat <- apply(kgrid, 1, function(kk) {
    exp(1i * pi * colSums(kk * t(xs)))
  })
  dmat1 <- t(t(mat) * kgrid[, 1] * 1i * pi)
  dmat2 <- t(t(mat) * kgrid[, 2] * 1i * pi)
  vals <- (mat %*% zs)[, 1]
  d1 <- (dmat1 %*% zs)[, 1]
  d2 <- (dmat2 %*% zs)[, 1]
  list(vals = vals, d1 = d1, d2 = d2)
}

circ_confine <- function(xs, vals, d1, d2, pow = 1) {
  xtx <- rowSums(xs^2)
  gfunc <- (xtx^pow - 1)^2
  g1 <- 4 *(xtx - 1) * pow * xtx^(pow-1) * xs[, 1]
  g2 <- 4 *(xtx - 1) * pow * xtx^(pow-1) * xs[, 2]
  gfunc[xtx > 1] <- 0;
  g1[xtx > 1] <- 0; g2[xtx > 1] <- 0  
  new_vals <- vals * gfunc
  new_d1 <- vals * g1 + gfunc * d1
  new_d2 <- vals * g2 + gfunc * d2
  list(vals = new_vals, d1 = new_d1, d2 = new_d2)
}

ctp <- function(vals, d1, d2) {
  contour(xseq, xseq, matrix(Re(vals), nrow=length(xseq) ), lwd = 2)
  contour(xseq, xseq, matrix(Re(d1), nrow=length(xseq)) , add=TRUE, col="red")
  contour(xseq, xseq, matrix(Re(d2), nrow=length(xseq) ), add=TRUE, col="blue")  
}

xseq <- seq(-1, 1, by = 0.01)
xs <- cbind(rep(xseq, length(xseq)), rep(xseq, each = length(xseq)))
lambda <- 8
k.basis <- 8
res_grf <- grf_grid(xs, lambda, k.basis)
res_grf %$% ctp(vals, d1, d2)

res_c <- res_grf %$% circ_confine(xs, vals, d1, d2, pow = 5)
res_c %$% ctp(vals, d1, d2)

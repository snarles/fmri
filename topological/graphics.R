library(lineId)
library(magrittr)

grf_zs <- function(lambda = 2, k.basis = 5) {
  kseq <- (-k.basis):k.basis
  kgrid <- cbind(rep(kseq, each=(2*k.basis + 1)), rep(kseq, (2*k.basis + 1)))
  klen <- dim(kgrid)[1]
  zs <- rnorm(klen) + 1i * rnorm(klen)
  zs <- zs * exp(-lambda/4 * rowSums(kgrid^2))
  zs  
}

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
  list(vals = vals, d1 = d1, d2 = d2, zs = zs)
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

cc_query_ <- function(zs, k.basis, pow) {
  kseq <- (-k.basis):k.basis
  kgrid <- cbind(rep(kseq, each=(2*k.basis + 1)), rep(kseq, (2*k.basis + 1)))
  ansf <- function(x) {
    ## rf vals and ds
    ef_vals <- exp(1i * pi * kgrid %*% x)[, 1]
    rf_val <- sum(zs * ef_vals)
    rf_d1 <- sum(1i*pi*kgrid[,1] * zs * ef_vals)
    rf_d2 <- sum(1i*pi*kgrid[,2] * zs * ef_vals)
    ## confine vals and ds
    xtx <- sum(x^2)
    gfunc <- (xtx^pow - 1)^2
    g1 <- 4 *(xtx - 1) * pow * xtx^(pow-1) * x[1]
    g2 <- 4 *(xtx - 1) * pow * xtx^(pow-1) * x[2]
    val <- gfunc * rf_val
    d1 <- gfunc * rf_d1 + g1 * rf_val
    d2 <- gfunc * rf_d2 + g2 * rf_val
    c(val = val, d1 = d1, d2 = d2)
  }
  ansf
}

## GENERATE X GRID
xseq <- seq(-1, 1, by = 0.01)
xs <- cbind(rep(xseq, length(xseq)), rep(xseq, each = length(xseq)))

## GENERATE RF
lambda <- 8
k.basis <- 8
res_grf <- grf_grid(xs, lambda, k.basis)
res_grf %$% ctp(vals, d1, d2)

## CONFINE RF TO CIRCLE
pow <- 5
res_c <- res_grf %$% circ_confine(xs, vals, d1, d2, pow)
res_c %$% ctp(vals, d1, d2)

## CHECK QUERY FUNC
query <- cc_query_(res_grf$zs, k.basis, pow)
ind <- sample(which(rowSums(xs^2) < 1), 1)
xs[ind, ]
res_grf %$% c(vals[ind], d1[ind], d2[ind])
res_c %$% c(vals[ind], d1[ind], d2[ind])
query(xs[ind, ])

## APPLY measure-preserving FLOW TO SPIRAL
ts <- seq(0, 1, by = 0.001) * 2 * pi
pts <- ((ts/pi)-1) * cbind(cos(20 * ts), sin(20 * ts))
par(bg = "grey")
plot(pts, xlim = c(-1,1), ylim = c(-1, 1), type = "l", lwd = 2)
eps <- 0.01
nits <- 50
for (i in 1:nits) {
  ds <- t(apply(pts, 1, query))
  pts <- pts + eps * Re(cbind(ds[,3], -ds[, 2]))
  if (i %% 10==0) lines(pts, col = rainbow(nits)[i], lwd = 2)
}

## APPLY non-measure-preserving FLOW TO SPIRAL
ts <- seq(0, 1, by = 0.001) * 2 * pi
pts <- ((ts/pi)-1) * cbind(cos(20 * ts), sin(20 * ts))
par(bg = "grey")
plot(pts, xlim = c(-1,1), ylim = c(-1, 1), type = "l", lwd = 2)
eps <- 0.01
nits <- 50
for (i in 1:nits) {
  ds <- t(apply(pts, 1, query))
  pts <- pts + eps * Re(cbind(ds[,2], ds[, 3]))
  if (i %% 10==0) lines(pts, col = rainbow(nits)[i], lwd = 2)
}

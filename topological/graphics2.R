source("topological//graphics_source.R")
library(akima)

xseq <- seq(-1, 1, by = 0.01)
xs <- cbind(rep(xseq, length(xseq)), rep(xseq, each = length(xseq)))
lambda <- 2
k.basis <- 8
res_c <- grf_grid(xs, lambda, k.basis)
xs0 <- xs;xs0 <- xs0[rowSums(xs^2) < 1, ]
fs <- Re(res_c$vals); fs <- fs[rowSums(xs^2) < 1]
eps1 <- matrix(Re(grf_grid(xs, lambda, k.basis)$vals), nrow = length(xseq))

xiseq <- seq(-1, 1, by = 0.01)

## Initial distortion

eps <- 0.01; nits <- 10; mp <- TRUE

xs1 <- list(xs0)
for (i in 1:2) {
  zs <- grf_zs(lambda, k.basis)
  query <- cc_query_(zs, k.basis, pow)
  xs1 <- deform_map(xs1, query, eps, nits, mp)
}
xs1 <- xs1[[1]]

## Redistortion

eps <- 0.01; nits <- 1; mp <- TRUE

xs2 <- list(xs1)
for (i in 1:2) {
  zs <- grf_zs(lambda, k.basis)
  query <- cc_query_(zs, k.basis, pow)
  xs2 <- deform_map(xs2, query, eps, nits, mp)
}
xs2 <- xs2[[1]]


is.lim <- function(a) sum(is.na(a))==0 & (sum(a^2) < 1)
inds <- apply(xs2, 1, is.lim)
fld <- interp(x = xs2[inds, 1], y = xs2[inds, 2], z = fs[inds],
              xo = xiseq, yo = xiseq)
z <- fld$z
z[is.na(z)] <- 0
xr <- matrix(fld$x[row(z)], nrow=nrow(z))
yr <- matrix(fld$y[col(z)], nrow=nrow(z))
## standardize z inside circle
zmax <- max(z[(xr^2 + yr^2) < 1])
zmin <- min(z[(xr^2 + yr^2) < 1])
z <- (z - zmin)/(zmax - zmin)
## set z outside to 0
z[(xr^2 + yr^2 >= 1)] <- 0

## Plot it

par(bg = "black")
ts <- seq(0, 1, by = 0.001) * 2 * pi
circ0 <- cbind(cos(ts), sin(ts))
plot(circ0, xlim = c(-1,1), ylim = c(-1, 1), type = "l", 
     lwd = 3, col = "white", axes = FALSE)
polygon(circ0, col="white")
.filled.contour(x = fld$x,
                y = fld$y,
                z = z,
                levels = seq(0, 1, length.out = 40),
                col = gray(0:39/40)
)
polygon(circ0, col=NA, border = "white", lwd = 5)





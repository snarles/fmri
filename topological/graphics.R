source("topological//graphics_source.R")
library(akima)

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


## Poly circle deformations

polys <- polycirc(10, radmax = 0.3, innmax = 0.7)[[1]]
par(bg = "grey");plot(NA, xlim = c(-1,1), ylim = c(-1, 1))
plot_polys(polys, col = "black")
nits <- 10
for (i in 1:nits) {  
  polys <- rand_deform(polys, lambda, k.basis, eps=0.01)
  plot_polys(polys, col = rainbow(nits)[i])  
}

def_plot <- function(bgc = "white", circ = "black", fc = "grey", ...) {
  ts <- seq(0, 1, by = 0.001) * 2 * pi
  circ0 <- cbind(cos(ts), sin(ts))
  par(bg = bgc)
  plot(circ0, xlim = c(-1,1), ylim = c(-1, 1), type = "l", lwd = 3, col = circ, ...)
  polygon(circ0, col=fc)
}

## Poly circle, two shade
res <- polycirc(20, radmax = 0.4, innmax = 0.8, rad_dec=0.2,
                ina = 3, inb = 1)
polys1 <- res[[1]]; polys2 <- res[[2]]
eps <- 0.01; nits <- 10; mp <- TRUE
lambda <- 8; k.basis <- 8

for (i in 1:2) {
  zs <- grf_zs(lambda, k.basis)
  query <- cc_query_(zs, k.basis, pow)
  polys1 <- deform_map(polys1, query, eps, nits, mp)
  polys2 <- deform_map(polys2, query, eps, nits, mp)  
}
def_plot(bgc = "black", circ = "blue", fc = gray(0.1))
plot_polys(polys1, col = gray(0.2), border = NA)
plot_polys(polys2, col = gray(0.4), border = NA)

## plot GRF, distorted

xseq <- seq(-1, 1, by = 0.01)
xs <- cbind(rep(xseq, length(xseq)), rep(xseq, each = length(xseq)))
lambda <- 8
k.basis <- 8
res_grf <- grf_grid(xs, lambda, k.basis)
pow <- 5
res_c <- res_grf %$% circ_confine(xs, vals, d1, d2, pow)

xs1 <- list(xs)
for (i in 1:2) {
  zs <- grf_zs(lambda, k.basis)
  query <- cc_query_(zs, k.basis, pow)
  xs1 <- deform_map(xs1, query, eps, nits, mp)
}
xs1 <- xs1[[1]]

is.lim <- function(a) sum(is.na(a))==0 & (sum(a^2) < 1)
inds <- apply(xs1, 1, is.lim)
fld <- interp(x = xs1[inds, 1], y = xs1[inds, 2], z = Re(res_c$vals)[inds])
fld$z[is.na(fld$z)] <- 0

par(bg = "black")
def_plot(bgc="black", axes = FALSE)
.filled.contour(x = fld$x,
                y = fld$y,
                z = fld$z,
                levels = seq(min(fld$z), max(fld$z), length.out = 20),
                col = gray(0:19/20)
)



##

pts <- xs1[inds, ][1:10, ]
pp <- list(pts)
zs <- grf_zs(lambda, k.basis)
zs[1:10]
deform_map(pp, zs, eps, 10 * nits, mp)

##


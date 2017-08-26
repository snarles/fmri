library(tripack)
library(pracma)

fcirc <- function(x, y, rad, ...) {
  ts <- 0:20/20
  polygon(x + rad * cos(2*pi*ts), y + rad * sin(2*pi*ts), ...)
}

par(bg = "white")
K <- 7; r <- 20; n <- K * r; sigma <- 0.2
xs <- randn(K, 2)
zs <- rep(1:K, each = r)
ys <- xs[zs, ] + sigma * randn(n, 2)
plot(NA, NA , axes = FALSE, cex = 3, xlim = c(-3, 3), ylim = c(-3, 3), ann = FALSE, pch = "+")
cols <- rainbow(K)
for (i in 1:n) {
  fcirc(ys[i, 1], ys[i, 2], 0.07, col = "grey", border = cols[zs[i]], lwd = 4)
}
for (i in 1:K) {
  fcirc(xs[i, 1], xs[i, 2], 0.25, col = "black", border = cols[i], lwd = 4)
  points(xs[i, , drop = FALSE], pch = paste(i), col = cols[i], cex = 2)
}

vm <- voronoi.mosaic(xs[, 1], xs[,2]) 
plot(vm, add = TRUE)

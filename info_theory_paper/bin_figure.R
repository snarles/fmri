####
##  Interpreting mutual information in terms of binning
####

library(MASS)

## GENERATE DATA

rho <- 0.995
n <- 10000
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
YX <- mvrnorm(n, c(0, 0), Sigma)
#YX[, 2] <- sign(YX[, 2]) * sqrt(abs(YX[, 2]))
#YX[, 1] <- sign(YX[, 2]) * sqrt(abs(YX[, 1]))
o <- order(YX[, 2])
YX <- YX[o, ]
ts <- (1:n)/n

K <- 10
#breaks <- min(YX[, 1])
breaks <- -1.5
ds <- list()
for (i in 1:K) {
  sub <- (ts <= i/K) & (ts > (i-1)/K)
  breaks <- c(breaks, max(YX[sub, 2]))
  Ysamp <- YX[sub, 1]
  ds[[i]] <- density(Ysamp)
}


## PLOTTING

plot(YX, xlim = c(-3, 3), ylim = c(-2, 2), ann = FALSE,
     axes = FALSE, cex.lab = 3, col = rgb(0,0,0, 0.1))
title(xlab = "Y", mgp = c(2, 1, 0.1), cex.lab = 3)
title(ylab = "X", mgp = c(0.5, 1, 0.1), cex.lab = 3)
axis(1, at = c(-3, 3), labels = FALSE, lwd = 2)
axis(2, at = c(-3, 3), labels = FALSE, lwd = 2)
for (i in 2:K) abline(breaks[i], 0, lwd = 2)

plot(NA, NA, xlim = c(-3, 3), ylim = c(-2, 2), ann = FALSE,
     axes = FALSE, cex.lab = 3, col = rgb(0,0,0, 0.1))
title(xlab = "Y", mgp = c(2, 1, 0.1), cex.lab = 3)
title(ylab = "X", mgp = c(0.5, 1, 0.1), cex.lab = 3)
axis(1, at = c(-3, 3), labels = FALSE, lwd = 2)
axis(2, at = c(-3, 3), labels = FALSE, lwd = 2)
for (i in 1:K) {
  polygon(ds[[i]]$x, 0.3 * ds[[i]]$y + breaks[i], col = rgb(0, 0, 0, 0.5), lwd = 2)
}
for (i in 1:K) abline(breaks[i], 0, lwd = 2)

-1/2 * log(1 - rho^2)

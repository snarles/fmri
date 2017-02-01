source("approximation/gaussian_lc_source.R")
library(mvtnorm)

rho <- 0.3
Sigma2 <- 1/(1-rho^2)
Sigma <- matrix(c(1,rho, rho, 1), 2, 2)


udelt <- 0.01
us <- seq(0, 1, udelt)

delta <- 0.01
xmax <- 4
dmax <- 7
xgrid <- seq(0, xmax, delta)
#xygrid <- cbind(rep(xgrid, each = length(xgrid)), rep(xgrid, length(xgrid)))
dgrid <- seq(0, dmax, delta)
xdgrid <- cbind(rep(xgrid, each = length(dgrid)), rep(dgrid, length(xgrid)))
Phi_diffs <- 1 - (pnorm(rho*xdgrid[, 1] + xdgrid[, 2]) - 
                    pnorm(rho*xdgrid[, 1] - xdgrid[, 2]))
Phi_diffs2 <- (1-(pnorm(rho*xdgrid[, 1] + xdgrid[, 2], mean = rho * xdgrid[, 1],
                    sd = sqrt(1-rho^2))-
                   pnorm(rho*xdgrid[, 1] - xdgrid[, 2], mean = rho * xdgrid[, 1],
                         sd = sqrt(1-rho^2)))) * dnorm(xdgrid[, 1])
#tab_all <- data.frame(x = xdgrid[, 1], d = xdgrid[, 2], Phi_diffs, Phi_diffs2);View(tab_all)
xP2 <- cbind(xdgrid[, 1], Phi_diffs2)

plot(NA, NA, xlim = range(xgrid), ylim = range(dgrid), xlab = "x", ylab = "d")
for (u in rev(seq(0, 1, 0.01))) {
  points(xdgrid[Phi_diffs <= u & Phi_diffs > (u - 0.01), ], pch = ".", col = grey(u))
}

Kfunc <- sapply(us, function(u) {
  xp <- xP2[Phi_diffs <= u, ]
  mx <- match(xgrid, xp[, 1])
  p2 <- xp[mx, 2]
  2 * sum(p2 * delta)
})
Kfunc[us==0] <- 0

plot(us, Kfunc, type = "l", ylim = c(0, 1))


for (k in 2:5) {
  print(c(emp=mc(matrix(Sigma2,1,1), k, mc.reps = 1e4),
          the=(k-1) * sum(Kfunc * us^(k-2)) * udelt))
}

dKfunc <- (Kfunc[-1] - Kfunc[-length(Kfunc)])/udelt
plot(us[-1] + udelt/2, dKfunc, type = "l")

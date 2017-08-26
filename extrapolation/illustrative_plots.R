

## plot mc errors
des <- readRDS("extrapolation/bc_rho0.7.rds")
plot(NA, NA, xlim = c(1, 11), ylim = c(0, 1), 
     xlab = "k", ylab = "Error", main = "rho = 0.7")
scale_fac <- 1/20
for (k in 2:10) {
  dd <- des[[paste(k)]]
  lines(scale_fac * dd$y + k, dd$x)
  lines(-scale_fac * dd$y + k, dd$x)
  points(k, sum(dd$x * dd$y)/sum(dd$y))
}

kures <- readRDS("extrapolation/ku_rho0.7.rds")
lines(kures$tabl[, "k"], kures$tabl[, "the.l"], col = "red", lwd = 2)

## plot of bivariate density
rho <- 0.7
Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
S.5 <- pracma::sqrtm(Sigma)$B
thetas <- seq(0, 2 * pi, pi/100)
circ <- cbind(cos(thetas), sin(thetas))
nqs <- qnorm(seq(0, 0.5, 0.05))
plot(0, 0, xlim = c(-2.5,2.5), asp = 1, col = "white",
     xlab = "x", ylab = "y", main = "rho = 0.7")
for (i in rev(1:9)) {
  polygon(circ %*% S.5 * qnorm(0.5 + i/20), 
          border = "black",
          col = grey(i/10))
}
abline(0, 1/rho, lty = 2, col = "red", lwd = 2)
legend(0.37, -2.4, col = "red", lty = 2, lwd = 2, legend = "E[X|Y=y]       ")

## next plot adds onto ^^
polygon(c(-9, 9, 9, -9), c(-9, -9, 9, 9), col = rgb(1,1,1,0.7))

delta <- 0.01
xmax <- 3.5
dmax <- 10
xgrid <- seq(0, xmax, delta)
dgrid <- seq(0, dmax, delta)
xdgrid <- cbind(rep(xgrid, each = length(dgrid)), rep(dgrid, length(xgrid)))
Phi_diffs <- 1 - (pnorm(xdgrid[, 1]/rho + xdgrid[, 2]) - 
                    pnorm(xdgrid[, 1]/rho - xdgrid[, 2]))
tab_all <- data.frame(x = xdgrid[, 1], d = xdgrid[, 2], Phi_diffs, 
                      lb = xdgrid[, 1]/rho - xdgrid[, 2],
                      ub = xdgrid[, 1]/rho + xdgrid[, 2])
us <- seq(0.1, 0.9, by = 0.1)
for (u in us) {
  stab <- subset(tab_all, Phi_diffs <= u)
  mx <- match(xgrid, stab$x)
  stab2 <- stab[mx, ]
  lines(stab2$x, stab2$lb)
  lines(stab2$x, stab2$ub)
  lines(-stab2$x, -stab2$lb)
  lines(-stab2$x, -stab2$ub)
  text(2, mean(stab2$lb[abs(stab2$x - 2) < 0.1]) + 0.11, 
       paste("U > ", u))
  text(-2, -mean(stab2$lb[abs(stab2$x - 2) < 0.1]) - 0.13, 
       paste("U > ", u))
  
}

## Ku function

plot(kures$us, kures$Kfunc, type = "l", ylim = c(0, 1),
     xlab = "u", ylab = "K(u)", main = "rho = 0.7")


## plot Kus and errs
dsets <- list()
rhos <- seq(0.3, 0.9, 0.2)
for (rho in rhos) {
  dsets[[paste(rho)]] <- readRDS(paste0("extrapolation/ku_rho", rho, ".rds"))
}

cols <- c("black", "green", "red", "blue"); names(cols) <- names(dsets)

plot(NA, NA, xlim = c(0,1), ylim = c(0,1), xlab = "u", ylab = "K(u)")
for (rho in rhos) {
  lines(dsets[[paste(rho)]]$us, dsets[[paste(rho)]]$Kfunc,
        lwd = 2, col = cols[paste(rho)])
}
legend(0, 1, col = cols, legend = paste("rho = ", rhos, "  "), lwd = 2)

plot(NA, NA, xlim = c(2,10), ylim = c(0,1), xlab = "k", ylab = "AvRisk")
for (rho in rhos) {
  lines(dsets[[paste(rho)]]$tabl[, "k"], dsets[[paste(rho)]]$tabl[, "the.l"],
        lwd = 2, col = cols[paste(rho)], type = "o")
}
legend(6, 0.3, col = cols, legend = paste("rho = ", rhos, "  "), lwd = 2, pch = "o")

## get polynomial approximation rates
source("extrapolation/polynomial_approximation.R")
dmax <- 10
tabl <- matrix(0, dmax, length(rhos))
colnames(tabl) <- paste(rhos)

for (rho in rhos) {
  ds <- dsets[[paste(rho)]]
  for (d in 1:dmax) {
    res <- uniform_poly(d, ds$us, ds$Kfunc)
    tabl[d, paste(rho)] <- res$error
  }
}

tabl
matplot(log(tabl), type = "o", ylab = "log(error)",
        xlab = "deg", lwd = 2, main = "Approximation Errors")
legend(6, -2.5, legend = paste("rho =", rhos, "  "), col = 1:4, lty = 1:4, lwd = 2, pch = paste(1:4))

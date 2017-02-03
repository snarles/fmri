

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

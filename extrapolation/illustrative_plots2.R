

## demonstrate how ku-produces errors
layout(matrix(1:2, 1, 2))

des <- readRDS("extrapolation/bc_rho0.7.rds")

kures <- readRDS("extrapolation/ku_rho0.7.rds")

for (kk in 2:10) {
  png(paste0("extrapolation/", "rho_0_7_fmla", kk, ".png"),
      width = 700, height = 450)
  layout(matrix(1:2, 1, 2))
  plot(NA, NA, xlim = c(1, 11), ylim = c(0, 1), 
       xlab = "k", ylab = "Error")
  scale_fac <- 1/20
  for (k in 2:10) {
    dd <- des[[paste(k)]]
    lines(scale_fac * dd$y + k, dd$x)
    lines(-scale_fac * dd$y + k, dd$x)
    points(k, sum(dd$x * dd$y)/sum(dd$y))
    if (k == kk) {
      points(k, sum(dd$x * dd$y)/sum(dd$y), 
             col = "red", cex = 3, pch = "+")
    }
  }
  plot(kures$us, kures$Kfunc, type = "l", ylim = c(0, 1),
       ylab = "D(u)", xlab = "u")
  bquants <- qbeta(seq(0.025, 0.975, 0.05), shape1 = kk-1,
                   shape2 = 1)
  quantpts <- approx(kures$us, kures$Kfunc, bquants)
  points(quantpts)
  polygon(c(kures$us, 1, 0), 
          c(0.1 * (kk-1) * kures$us^(kk-2),
            0, 0),
          border = "red", col = rgb(1, 0, 0, 0.5))
  dev.off()
}
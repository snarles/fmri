library(pracma)

reso <- 1000
xs <- seq(0, 1, length.out = reso)
qx <- (0.01 * xs^2 + 0.5 * exp(xs^20) - 0.5)^(1/3)
qx <- qx/sum(qx) * reso
#plot(xs, qx, type = "l")

par(mar = c(5,4,2,0))
layout(matrix(1:2, 1, 2))

mat <- matrix(0, reso, reso)
idx <- 1 + row(mat) - col(mat) + reso * (row(mat) < col(mat))
mat <- matrix(qx[idx], reso, reso)
plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "y")
image(flipud(t(mat)), col = rev(grey.colors(n = 100, start = 0, end = 1)),
      add = TRUE)
for (i in 0:4) abline(i/5, 0, col = "red", lty = 2)

par(mar = c(5,2,2,2))
plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "")
for (i in 1:5) {
  pos <- (i-1) * 200 + 1
  qx2 <- c(qx[pos:reso], qx[-(pos:reso)])
  lines(xs, qx2/30 + pos/reso)
  abline(pos/reso, 0, col = "red", lty = 2)
}
text(0.5, 0.15, "p(x|y=0) = q(x)")
text(0.5, 0.35, "p(x|y=0.2)")
text(0.2, 0.55, "p(x|y=0.4)")
text(0.6, 0.75, "p(x|y=0.6)")
text(0.5, 0.95, "p(x|y=0.8)")


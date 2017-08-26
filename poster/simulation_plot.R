load('results2.RData')
mus <- apply(results, c(2, 3), mean)
sds <- apply(results, c(2, 3), sd)
n <- 1:15 * 5 + 5
mus

pdf("poster/simulation1.pdf", width = 7, height = 4)
plot(n, mus[1, ], ylim = c(0, 1), ylab = "err", type = "o", lwd = 2, cex.lab = 2, cex.axis = 2)
lines(n, mus[4, ], col = "red", type = "o", pch = "M", lwd = 2)
lines(n, mus[5, ], col = "blue", type = "o", pch = "E", lwd = 2)
for (i in 1:15) {
  lines(n[i] * c(1, 1), mus[1, i] + c(1, -1) * sds[1, i])  
  points(n[i] * c(1, 1), mus[1, i] + c(1, -1) * sds[1, i], pch = "-")  
  lines(n[i] * c(1, 1), mus[2, i] + c(1, -1) * sds[2, i], col = "red")  
  points(n[i] * c(1, 1), mus[2, i] + c(1, -1) * sds[2, i], pch = "-", col = "red")
  lines(n[i] * c(1, 1), mus[5, i] + c(1, -1) * sds[5, i], col = "blue")  
  points(n[i] * c(1, 1), mus[5, i] + c(1, -1) * sds[5, i], pch = "-", col = "blue")  
}
dev.off()
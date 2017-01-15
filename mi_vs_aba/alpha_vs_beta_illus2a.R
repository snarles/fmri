## get the curve

source("mi_vs_aba/alpha_vs_beta.R")

for (k in 2:5) {
  reso <- 1000
  
  # pts <- sapply(1:10000, function(i) {
  #   qs <- rbeta(reso, runif(1) * 10, runif(1) * 10)
  #   qs <- normalize_qs(qs)
  #   get_q_ints(qs, k)
  # })
  bts <- 1000 * (1:10000/10000)
  pts2 <- sapply(bts, function(x) {
    qs <- qs_par(x, k, reso)
    get_q_ints(qs, k)
  })
  #plot(bts, pts2[2, ], type = "l")
  #plot(sqrt(bts), pts2[2, ], type = "l")
  # plot(log(bts), pts2[2, ], type = "l")
  # plot(bts, pts2[3, ], type = "l")
  # plot(log(bts), pts2[3, ], type = "l")
  
  # plot(t(pts)[, -1], pch= ".")
  # lines(t(pts2)[, -1], col = "red")
  # 
  # 
  # 
  # 
  # plot(qs_par(d = 2, k = 3), type = "l")
  # plot(qs_par(d = 2, k = 4), type = "l")
  # plot(qs_par(d = 2, k = 4), type = "l")
  # plot(qs_par(d = 10, k = 4), type = "l")
  # 
  
  ppts <- t(pts2)[, -1]
  ppts <- ppts[complete.cases(ppts), ]
  ppts2 <- cbind(1-ppts[, 2], ppts[, 1])
  plot(1 - ppts[, 2], ppts[, 1], type = "l", ylim = c(0, 5), 
       xlim = c(0, 1), xlab = "IdRisk", ylab = "I", lwd = 2)
  polygon(rbind(ppts2, c(1 - 1/k, max(ppts[, 1], na.rm = TRUE)), c(1-1/k, 0)), col = "grey", border = NA)
  lines(ppts2[, 1], ppts2[, 2], lwd =2)
  #text(0.5, ppts[, 2][which(ppts[, 1] > 0.5)[1]] + 0.1, expression(C[k]))
  # lines(c(0, max(ppts[, 1])), c(1/k, 1/k))
  lala <- expression(g[3])
  lala[[1]][[3]] <- k
  title(lala, cex.main = 2)
  # text(3.5, 0.5, "S", cex = 3)
  
}

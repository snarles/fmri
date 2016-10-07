## get the curve

source("mi_vs_aba/alpha_vs_beta.R")

for (k in 2:5) {
  reso <- 1000
  
  # pts <- sapply(1:10000, function(i) {
  #   qs <- rbeta(reso, runif(1) * 10, runif(1) * 10)
  #   qs <- normalize_qs(qs)
  #   get_q_ints(qs, k)
  # })
  
  pts2 <- sapply(1000 * (1:10000/10000), function(x) {
    qs <- qs_par(x, k, reso)
    get_q_ints(qs, k)
  })
  
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
  plot(ppts[, 1], ppts[, 2], type = "l", ylim = c(0, 1), 
       xlim = c(0, 5), xlab = "I", ylab = "ABA", lwd = 2)
  polygon(rbind(ppts, c(max(ppts[, 1], na.rm = TRUE), 1/k), c(0, 1/k)), col = "grey", border = NA)
  lines(ppts[, 1], ppts[, 2], lwd =2)
  #text(0.5, ppts[, 2][which(ppts[, 1] > 0.5)[1]] + 0.1, expression(C[k]))
  # lines(c(0, max(ppts[, 1])), c(1/k, 1/k))
  title(paste("k =", k))
  # text(3.5, 0.5, "S", cex = 3)
  
}
